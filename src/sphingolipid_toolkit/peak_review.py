"""Read-only chromatogram and spectrum review for an annotated local feature."""

import ast
from dataclasses import dataclass
import gzip
import hashlib
import json
import math
from pathlib import Path
from queue import Empty, Queue
import textwrap
import threading
import tkinter as tk
from tkinter import ttk

import numpy as np
from pyteomics import mzml

from .gui_theme import (ACCENT as MATCH, ALTERNATE, AXIS as GREY, BACKGROUND, HEADER, INK,
                        MUTED as TEXT, PANEL as PAPER, RULE, SELECTION as AREA, UI_FONT, style_review_axis)
from .targeted_mzml_pipeline import _extract_scan_rt


POOL_LABELS = {"hg": "HG", "lcb": "LCB", "other": "Other", "fah": "FAH", "b": "B", "loss": "Loss"}
PEAK_FIELDS = {
    "source_path", "polarity", "mz_min", "mz_max", "rt_min", "rt_max",
    "apex_rt_min", "apex_native_scan_id", "native_scan_id", "apex_intensity",
    "lipid_name", "observed_mz", "total_score", "matched_fragment_count",
}


def _peaks(scan):
    mass = np.asarray(scan.get("m/z array", []), dtype=float)
    intensity = np.asarray(scan.get("intensity array", []), dtype=float)
    keep = np.isfinite(mass) & np.isfinite(intensity) & (intensity > 0)
    return mass[keep], intensity[keep]


def read_peak_scans(row, source):
    """Use the same mass bounds, summed centroids and polarity as the peak export."""
    start, end = float(row.rt_min), float(row.rt_max)
    margin = max((end - start) * 0.65, 0.055)
    lo, hi = float(row.mz_min), float(row.mz_max)
    epsilon = max(abs(hi) * 2e-14, 1e-12)
    mode = {"+": "positive scan", "-": "negative scan"}[row.polarity]
    apex_id = str(row.apex_native_scan_id)
    points = []
    with mzml.MzML(str(source), use_index=True) as reader:
        ids = tuple(reader.default_index)
        if apex_id not in reader.default_index or row.native_scan_id not in reader.default_index:
            raise ValueError("The recorded scan IDs are missing from this mzML file.")
        apex_index = ids.index(apex_id)
        for indices in (range(apex_index, -1, -1), range(apex_index + 1, len(ids))):
            for index in indices:
                scan = reader.get_by_id(ids[index])
                rt = _extract_scan_rt(scan)
                if rt is None:
                    continue
                if rt < start - margin or rt > end + margin:
                    break
                if scan.get("ms level") == 1 and mode in scan:
                    mass, intensity = _peaks(scan)
                    within = (mass >= lo - epsilon) & (mass <= hi + epsilon)
                    points.append((rt, float(intensity[within].sum()), ids[index]))
        tandem = reader.get_by_id(row.native_scan_id)
        tandem_rt = _extract_scan_rt(tandem)
        if tandem.get("ms level") != 2 or mode not in tandem or tandem_rt is None:
            raise ValueError("The recorded tandem scan does not match this annotation.")
        if not start <= tandem_rt <= end:
            raise ValueError("The tandem scan is outside the recorded peak boundaries.")
        spectrum = _peaks(tandem)
    apex = next((point for point in points if point[2] == apex_id), None)
    if apex is None or not np.isclose(apex[0], row.apex_rt_min, rtol=0, atol=1e-7) or not np.isclose(
            apex[1], row.apex_intensity, rtol=1e-6, atol=1e-7):
        raise ValueError("The measured MS1 apex does not match the result workbook.")
    trace = np.array([(rt, value) for rt, value, _ in sorted(points)])
    return trace, spectrum, tandem_rt


@dataclass(frozen=True)
class MatchedIon:
    mz: float
    intensity: float
    name: str
    pool: str = ""


def _annotation_evidence(row, evidence_dir):
    """Resolve the workbook's evidence token without attaching another row's data."""
    raw = row.get("ms2_evidence_json")
    if not isinstance(raw, str):
        return None
    try:
        data = json.loads(raw)
        if not isinstance(data, dict):
            return None
        if "external_evidence_sha256" in data:
            name = row.get("evidence_file")
            if evidence_dir is None or not isinstance(name, str) or not name:
                return None
            directory = Path(evidence_dir).resolve()
            sidecar = (directory / name).resolve()
            if not sidecar.is_relative_to(directory) or sidecar.suffix.lower() != ".gz":
                return None
            key = data["external_evidence_sha256"]
            with gzip.open(sidecar, "rt", encoding="utf-8") as handle:
                raw = json.load(handle)[key]
            if hashlib.sha256(raw.encode("utf-8")).hexdigest() != key:
                return None
            data = json.loads(raw)
        return data.get("fragments") if isinstance(data, dict) else None
    except (OSError, ValueError, TypeError, AttributeError, KeyError, EOFError):
        return None


def _ion_values(row, evidence_dir):
    fragments = _annotation_evidence(row, evidence_dir)
    if isinstance(fragments, list) and all(isinstance(fragment, dict) for fragment in fragments):
        for fragment in fragments:
            if fragment.get("observed_mz") is not None:
                pool = POOL_LABELS.get(str(fragment.get("scoring_pool")), "")
                if pool == "Other" and fragment.get("auxiliary_lcb_support") is True:
                    pool = "Other (aux. LCB)"
                yield fragment["observed_mz"], fragment.get("intensity"), fragment.get("name", ""), pool
        return
    # Older workbooks also carry the three paired observation columns.
    values = []
    for key in ("实际mz", "强度", "碎片"):
        raw = row.get(key)
        try:
            array = ast.literal_eval(raw) if isinstance(raw, str) else raw
        except (ValueError, SyntaxError):
            return
        if not isinstance(array, (tuple, list)):
            return
        values.append(array)
    if len({len(array) for array in values}) != 1:
        return
    for mz, intensity, name in zip(*values):
        yield mz, intensity, name, ""


def _matched_ions(row, evidence_dir=None):
    ions = []
    for mass, abundance, name, pool in _ion_values(row, evidence_dir):
        try:
            mz, intensity = float(mass), float(abundance)
        except (TypeError, ValueError):
            continue
        if math.isfinite(mz) and math.isfinite(intensity) and mz > 0 and intensity > 0:
            ions.append(MatchedIon(mz, intensity, str(name), pool))
    return sorted(ions, key=lambda ion: ion.mz)


def _label_spectrum(axis, matches, mass, relative, base, limit):
    """Place a few prominent m/z labels above their own sticks."""
    half_width = limit * .074
    placed = []
    for ion in sorted(matches, key=lambda item: -item.intensity):
        height = 100 * ion.intensity / base
        if placed and height < 8:
            break
        if any(abs(ion.mz - mz) < .036*limit for mz, _, _, _ in placed):
            continue
        for align, left, right in (("center", ion.mz-half_width, ion.mz+half_width),
                                   ("left", ion.mz, ion.mz+2*half_width),
                                   ("right", ion.mz-2*half_width, ion.mz)):
            if left < 0 or right > limit:
                continue
            if np.any((left < mass) & (mass < right) & (relative > height + 5)):
                continue
            if any(left < end and start < right and abs(height - y) < 10 for _, start, end, y in placed):
                continue
            axis.annotate(f"{ion.mz:.4f}", (ion.mz, height), xytext=(0, 3.5), textcoords="offset points",
                          ha=align, va="bottom", color=MATCH, fontsize=8, zorder=5)
            placed.append((ion.mz, left, right, height))
            break
        if len(placed) == 4:
            break


def peak_figure(row, source, trace, spectrum, tandem_rt, *, matches=None):
    from matplotlib.figure import Figure
    from matplotlib.font_manager import fontManager
    from matplotlib.lines import Line2D
    from matplotlib.patches import Rectangle
    from matplotlib.ticker import MaxNLocator

    matches = _matched_ions(row) if matches is None else matches
    chart_top = 2.59
    chart_height = 2.30
    table_top = chart_top + chart_height + .53
    rows_top = table_top + .64
    names = [textwrap.wrap(ion.name, width=47) or [""] for ion in matches]
    row_heights = [.242 + .15 * (len(lines)-1) for lines in names]
    height = rows_top + max(.264, sum(row_heights)) + .17
    fig = Figure(figsize=(10, height), dpi=100, facecolor=PAPER)
    installed_fonts = {font.name for font in fontManager.ttflist}
    families = [UI_FONT[0], "DejaVu Sans"] if UI_FONT[0] in installed_fonts else ["DejaVu Sans"]
    for family in ("Microsoft YaHei", "Noto Sans CJK SC", "PingFang SC", "SimHei", "SimSun"):
        if family in installed_fonts:
            families.append(family)
            break

    def y_fraction(inches):
        return 1 - inches / height

    def label(x, y, text, size, *, color=INK, weight="normal", align="left"):
        fig.text(x, y_fraction(y), text, fontsize=size, color=color, weight=weight,
                 ha=align, va="center", family=families)

    def block(x, y, w, h, color):
        fig.add_artist(Rectangle((x, y_fraction(y+h)), w, h/height, transform=fig.transFigure,
                                 facecolor=color, edgecolor="none", zorder=0))

    def rule(y):
        y = y_fraction(y)
        fig.add_artist(Line2D([.045, .955], [y, y], transform=fig.transFigure,
                              color=RULE, linewidth=.65))

    block(0, 0, 1, .47, BACKGROUND)
    label(.048, .235, "SphinGOlipID", 12, weight="bold")
    label(.952, .235, "Review Peak & Spectrum", 10, color=TEXT, align="right")
    rule(.47)
    label(.048, .73, str(row.lipid_name), 16, weight="bold")
    mode = "Positive" if row.polarity == "+" else "Negative"
    label(.952, .73, f"{row.get('adduct', '')}   ·   {mode} ion mode", 10, align="right")
    label(.048, .99, Path(source).name, 9, color=TEXT)

    block(.045, 1.15, .91, .65, HEADER)
    metrics = [
        ("Precursor m/z", f"{row.observed_mz:.5f}"),
        ("Peak apex (min)", f"{row.apex_rt_min:.3f}"),
        ("MS2 score", f"{row.total_score:.2f}"),
        ("Matched peaks", str(int(row.matched_fragment_count))),
    ]
    for x, (name, value) in zip((.059, .301, .548, .788), metrics):
        label(x, 1.33, name, 9, color=TEXT)
        label(x, 1.60, value, 13, weight="bold")
    pool_scores = []
    for key, name in POOL_LABELS.items():
        try:
            value = float(row.get(f"{key}_score", 0))
        except (TypeError, ValueError):
            continue
        if math.isfinite(value) and value > 0:
            pool_scores.append(f"{name} {value:.2f}")
    if pool_scores:
        label(.048, 1.95, "Pool scores   ·   " + "   ·   ".join(pool_scores), 9, color=TEXT)
    rule(2.09)
    label(.048, 2.25, "MS1", 11, weight="bold")
    label(.525, 2.25, "MS2", 11, weight="bold")
    label(.048, 2.43, f"m/z {row.mz_min:.4f}–{row.mz_max:.4f}  ·  RT {row.rt_min:.3f}–{row.rt_max:.3f} min",
          9, color=TEXT)
    label(.525, 2.43, f"RT {tandem_rt:.3f} min", 9, color=TEXT)

    chart_y = y_fraction(chart_top + chart_height)
    left = fig.add_axes([.104, chart_y, .347, chart_height/height], facecolor=PAPER)
    right = fig.add_axes([.565, chart_y, .387, chart_height/height], facecolor=PAPER)
    for ax in (left, right):
        style_review_axis(ax)
        ax.set_yticks([0, 50, 100])
        ax.set_ylabel("Relative intensity (%)", fontsize=9, color=TEXT)
    rt, intensity = trace.T
    # Normalize against the entire displayed MS1 window. The recorded local apex
    # may be below 100% if a stronger neighbouring peak is visible.
    scale = max(float(intensity.max()), 1e-20)
    y = intensity * 100 / scale
    inside = (row.rt_min < rt) & (rt < row.rt_max)
    integrated_rt = np.r_[row.rt_min, rt[inside], row.rt_max]
    endpoint_y = np.interp([row.rt_min, row.rt_max], rt, y)
    integrated_y = np.r_[endpoint_y[0], y[inside], endpoint_y[1]]
    left.fill_between(integrated_rt, 0, integrated_y, facecolor=AREA, edgecolor="none", zorder=1)
    left.plot(rt, y, color=INK, lw=1.35, zorder=3)
    left.plot(row.apex_rt_min, row.apex_intensity * 100 / scale, "o", color=INK, ms=4, zorder=4)
    left.set_xlim(min(float(rt.min()), row.rt_min) - .001, max(float(rt.max()), row.rt_max) + .001)
    left.set_ylim(0, 112)
    left.xaxis.set_major_locator(MaxNLocator(nbins=4))
    left.ticklabel_format(axis="x", useOffset=False)
    left.set_xlabel("RT (min)", fontsize=10, color=TEXT)

    mass, counts = spectrum
    base = float(counts.max()) if counts.size else 0
    right.set_ylim(0, 119)
    if base:
        limit = max(float(mass.max()), row.observed_mz, max((ion.mz for ion in matches), default=0)) * 1.035
        right.set_xlim(0, limit)
        relative = 100 * counts / base
        right.vlines(mass, 0, relative, color=GREY, linewidth=.55, zorder=2)
        if matches:
            observed = [ion.mz for ion in matches]
            intensities = [100 * ion.intensity / base for ion in matches]
            right.vlines(observed, 0, intensities, color=MATCH, linewidth=1.5, zorder=3)
            _label_spectrum(right, matches, mass, relative, base, limit)
    right.xaxis.set_major_locator(MaxNLocator(nbins=5))
    right.set_xlabel("m/z", fontsize=10, color=TEXT)

    rule(table_top)
    label(.048, table_top+.19, "Matched ions", 11, weight="bold")
    if len(matches) != int(row.matched_fragment_count):
        label(.952, table_top+.19, f"{len(matches)} of {int(row.matched_fragment_count)} available in the workbook",
              9, color=TEXT, align="right")
    has_pools = any(ion.pool for ion in matches)
    mass_x = .746 if has_pools else .663
    block(.045, table_top+.36, .91, .28, HEADER)
    heading = table_top + .50
    label(.058, heading, "Ion assignment", 9, weight="bold")
    if has_pools:
        label(.525, heading, "Scoring pool", 9, weight="bold")
    label(mass_x, heading, "Observed m/z", 9, weight="bold", align="right")
    label(.94, heading, "Relative intensity (%)", 9, weight="bold", align="right")
    ypos = rows_top
    for index, (ion, lines, row_height) in enumerate(zip(matches, names, row_heights)):
        if index % 2:
            block(.045, ypos, .91, row_height, ALTERNATE)
        center = ypos + row_height/2
        label(.058, center, "\n".join(lines), 9.6)
        if has_pools:
            label(.525, center, ion.pool, 9.2)
        label(mass_x, center, f"{ion.mz:.4f}", 9.6, align="right")
        if base:
            label(.94, center, f"{100 * ion.intensity / base:.2f}", 9.6, align="right")
        ypos += row_height
        rule(ypos)
    return fig


def _show_peak_figure(win, figure):
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

    viewport = tk.Canvas(win, background=PAPER, highlightthickness=0, yscrollincrement=24)
    scrollbar = ttk.Scrollbar(win, orient=tk.VERTICAL, command=viewport.yview)
    viewport.configure(yscrollcommand=scrollbar.set)
    scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
    viewport.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
    width_inches, height_inches = figure.get_size_inches()
    canvas = FigureCanvasTkAgg(figure, viewport)
    window = viewport.create_window(0, 0, window=canvas.get_tk_widget(), anchor="nw")

    def fit_width(width):
        if width <= 1:
            return
        dpi = width / width_inches
        height = math.ceil(height_inches * dpi)
        figure.set_dpi(dpi)
        figure.set_size_inches(width_inches, height / dpi, forward=False)
        viewport.itemconfigure(window, width=width, height=height)
        viewport.configure(scrollregion=(0, 0, width, height))
        canvas.draw_idle()

    fit_width(max(1, win.winfo_width()-scrollbar.winfo_reqwidth()))
    viewport.bind("<Configure>", lambda event: fit_width(event.width))
    # TkAgg adjusts its DPI on mapping, including on scaled Windows displays.
    # Fit the plots to the visible client width after that adjustment.
    canvas.get_tk_widget().bind("<Map>", lambda _event: fit_width(viewport.winfo_width()), add="+")
    win.bind("<MouseWheel>", lambda event: viewport.yview_scroll(-int(event.delta/120), "units"))
    win.bind("<Button-4>", lambda _event: viewport.yview_scroll(-1, "units"))
    win.bind("<Button-5>", lambda _event: viewport.yview_scroll(1, "units"))


def open_peak_review(parent, row, source, *, evidence_dir=None):
    """Load raw scans off the Tk thread and draw only while the review is open."""
    win = tk.Toplevel(parent)
    win.title("SphinGOlipID · Review Peak & Spectrum")
    win.configure(background=PAPER)
    width = min(1000, win.winfo_screenwidth() - 60)
    height = min(930, win.winfo_screenheight() - 100)
    win.geometry(f"{width}x{height}")
    status = ttk.Label(win, text="Reading the selected peak from mzML…", anchor="center")
    status.pack(expand=True, fill="both")
    completed = Queue()

    def load():
        try:
            completed.put((read_peak_scans(row, source), _matched_ions(row, evidence_dir)))
        except Exception as exc:
            completed.put(exc)

    def display():
        nonlocal timer
        try:
            data = completed.get_nowait()
        except Empty:
            timer = win.after(75, display)
            return
        timer = None
        try:
            if isinstance(data, Exception):
                raise data
            scans, matches = data
            figure = peak_figure(row, source, *scans, matches=matches)
            _show_peak_figure(win, figure)
            status.destroy()
        except Exception as exc:
            if status.winfo_exists():
                status.configure(text=f"Could not display this feature.\n{exc}", style="Warning.TLabel")

    def cancel_poll(event):
        if event.widget is win and timer is not None:
            win.after_cancel(timer)

    timer = win.after(75, display)
    win.bind("<Destroy>", cancel_poll)
    threading.Thread(target=load, daemon=True).start()
    return win
