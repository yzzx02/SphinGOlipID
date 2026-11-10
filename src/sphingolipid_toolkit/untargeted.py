"""XCMS local MS1 features linked to the unchanged current MS2 identifier."""
from __future__ import annotations

from dataclasses import asdict, dataclass
from importlib.resources import files
import json
import logging
import math
import os
from pathlib import Path
import shutil
import subprocess
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

from .engine import IdentificationConfig, SphingolipidIdentifier, compatible_frame, library_path
from .reference_pipeline import iter_prepared_mzml
from .result_writer import write_identification_workbook


@dataclass(frozen=True)
class UntargetedConfig:
    mzml_files: tuple[Path, ...]
    output_file: Path
    polarity: str = "+"
    precursor_ppm: float = 10.0
    fragment_ppm: float | None = 20.0
    fragment_da: float | None = None
    min_total_score: float = 50.0
    top_n: int = 1
    peakwidth_min: float = 5.0
    peakwidth_max: float = 60.0
    noise: float = 1000.0
    snthresh: float = 3.0
    save_intermediate: bool = False
    rscript_path: Path | None = None

    def __post_init__(self):
        object.__setattr__(self, "mzml_files", tuple(Path(p).resolve() for p in self.mzml_files))
        object.__setattr__(self, "output_file", Path(self.output_file).resolve())
        if self.rscript_path is not None:
            object.__setattr__(self, "rscript_path", Path(self.rscript_path))
        if not self.mzml_files or len(set(self.mzml_files)) != len(self.mzml_files):
            raise ValueError("请选择至少一个 mzML 文件，且不要重复选择。")
        if any(p.suffix.lower() != ".mzml" for p in self.mzml_files):
            raise ValueError("输入必须是 mzML 文件。")
        if self.output_file.suffix.lower() != ".xlsx":
            raise ValueError("输出文件须为 .xlsx。")
        if self.output_file in self.mzml_files:
            raise ValueError("输出文件不能覆盖原始数据。")
        if self.polarity not in ("+", "-"):
            raise ValueError("请选择正模式或负模式。")
        for key in ("peakwidth_min", "peakwidth_max", "snthresh"):
            value = getattr(self, key)
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"{key} 须为大于零的有限数值。")
        if self.peakwidth_min >= self.peakwidth_max:
            raise ValueError("最小峰宽须小于最大峰宽，单位为秒。")
        if not math.isfinite(self.noise) or self.noise < 0:
            raise ValueError("噪声阈值须为非负有限数值。")
        self.identification_config()

    def identification_config(self):
        return IdentificationConfig(precursor_ppm=self.precursor_ppm,
            fragment_ppm=self.fragment_ppm, fragment_da=self.fragment_da,
            min_total_score=self.min_total_score, top_n=self.top_n)


def find_rscript(selected: Path | None = None) -> Path:
    if selected is not None:
        if not selected.is_file():
            raise FileNotFoundError(f"找不到 Rscript：{selected}")
        return selected.resolve()
    executable = shutil.which("Rscript")
    if executable:
        return Path(executable).resolve()
    if os.name == "nt":
        root = Path(os.environ.get("ProgramFiles", r"C:\Program Files")) / "R"
        candidates = list(root.glob("R-*/bin/Rscript.exe")) + list(root.glob("R-*/bin/x64/Rscript.exe"))
        if candidates:
            return max(candidates, key=lambda path: path.stat().st_mtime_ns)
    raise RuntimeError("未找到 Rscript。请安装 R、xcms、MSnbase、mzR、jsonlite 和 BiocParallel，并将 Rscript 加入 PATH。")


def run_xcms(path: Path, config: UntargetedConfig, work: Path, logger: logging.Logger):
    """Run one file serially; never align, fill or infer features from MS2."""
    work.mkdir(parents=True, exist_ok=True)
    request = work / "request.json"
    output, metadata = work / "chromatographic_peaks.csv", work / "xcms_versions.json"
    values = {key: getattr(config, key) for key in (
        "precursor_ppm", "peakwidth_min", "peakwidth_max", "noise", "snthresh", "polarity")}
    values.update(input=str(path), output=str(output), metadata=str(metadata))
    payload = json.dumps(values, ensure_ascii=True)
    request.write_text(payload, encoding="utf-8")
    script = Path(str(files("sphingolipid_toolkit").joinpath("r", "xcms_centwave.R")))
    # Rscript on Windows can misinterpret a Unicode configuration *filename*
    # before jsonlite opens it. Send ASCII-escaped JSON through standard input;
    # the decoded source/output paths retain their actual Unicode identity.
    command = [str(find_rscript(config.rscript_path)), "--vanilla", str(script), "--stdin"]
    logger.info("Detecting MS1 peaks with XCMS: %s", path.name)
    options = {"creationflags": subprocess.CREATE_NO_WINDOW} if os.name == "nt" else {}
    tail = []
    with subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          text=True, encoding="utf-8", errors="replace", **options) as process:
        try:
            process.stdin.write(payload)
            process.stdin.close()
        except BrokenPipeError:
            pass  # Read the actual R/package error below before reporting it.
        for line in process.stdout:
            text = line.rstrip()
            if text:
                tail.append(text)
                tail = tail[-25:]
                logger.info("XCMS: %s", text)
        status = process.wait()
    if status:
        raise RuntimeError("XCMS 峰提取失败：\n" + "\n".join(tail))
    if not output.is_file() or not metadata.is_file():
        raise RuntimeError("XCMS 未生成完整峰表和运行信息。")
    frame = pd.read_csv(output, encoding="utf-8-sig")
    return normalize_xcms_features(frame, path, config.polarity), json.loads(metadata.read_text(encoding="utf-8"))


FEATURE_COLUMNS = ["source_file", "source_path", "polarity", "xcms_peak_id", "feature_mz",
    "mz_min", "mz_max", "xcms_rt_min", "rt_min", "rt_max", "peak_area", "xcms_peak_height", "snr"]


def normalize_xcms_features(frame: pd.DataFrame, path: Path, polarity: str) -> pd.DataFrame:
    mapping = {"mz": "feature_mz", "mzmin": "mz_min", "mzmax": "mz_max", "rt": "xcms_rt_min",
               "rtmin": "rt_min", "rtmax": "rt_max", "into": "peak_area", "maxo": "xcms_peak_height", "sn": "snr"}
    if not {"xcms_peak_id", *mapping}.issubset(frame):
        raise ValueError("XCMS 峰表缺少必要的本地峰编号或测量字段。")
    frame = frame.copy().rename(columns=mapping)
    if frame.xcms_peak_id.isna().any() or frame.xcms_peak_id.duplicated().any():
        raise ValueError("XCMS 本地峰编号缺失或重复。")
    for name in mapping.values():
        frame[name] = pd.to_numeric(frame[name], errors="raise")
    required = frame[["feature_mz", "mz_min", "mz_max", "xcms_rt_min", "rt_min", "rt_max", "peak_area", "xcms_peak_height"]]
    if not np.isfinite(required.to_numpy()).all():
        raise ValueError("XCMS 返回非有限的峰测量值。")
    if ((frame.mz_min > frame.feature_mz) | (frame.feature_mz > frame.mz_max) |
        (frame.rt_min > frame.xcms_rt_min) | (frame.xcms_rt_min > frame.rt_max) |
        (frame.mz_min <= 0) | (frame.rt_min < 0) | (frame.peak_area < 0) | (frame.xcms_peak_height < 0)).any():
        raise ValueError("XCMS 返回无效的峰边界或强度。")
    for name in ("xcms_rt_min", "rt_min", "rt_max"):
        frame[name] /= 60.0
    frame["source_file"], frame["source_path"], frame["polarity"] = path.name, str(path.resolve()), polarity
    return frame[FEATURE_COLUMNS].reset_index(drop=True)


def measure_local_apices(features: pd.DataFrame, path: Path, polarity: str) -> pd.DataFrame:
    """Measure EIC maxima inside the actual XCMS mass/time bounds, streaming MS1."""
    from pyteomics import mzml
    from .targeted_mzml_pipeline import _extract_scan_rt
    result = features.copy()
    result["apex_rt_min"] = np.nan
    result["apex_mz"] = np.nan
    result["apex_intensity"] = np.nan
    result["apex_native_scan_id"] = ""
    if result.empty:
        return result
    starts = result.rt_min.to_numpy()
    ends = result.rt_max.to_numpy()
    mass_min, mass_max = result.mz_min.to_numpy(), result.mz_max.to_numpy()
    maximum = np.zeros(len(result))
    with mzml.MzML(str(path), decode_binary=True) as reader:
        for raw in reader:
            if raw.get("ms level") != 1 or ("positive scan" if polarity == "+" else "negative scan") not in raw:
                continue
            rt = _extract_scan_rt(raw)
            if rt is None:
                continue
            active = np.flatnonzero((starts <= rt) & (rt <= ends))
            if not len(active):
                continue
            masses = np.asarray(raw.get("m/z array", []), dtype=float)
            intensities = np.asarray(raw.get("intensity array", []), dtype=float)
            valid = np.isfinite(masses) & np.isfinite(intensities) & (intensities > 0)
            masses, intensities = masses[valid], intensities[valid]
            order = np.argsort(masses, kind="stable")
            masses, intensities = masses[order], intensities[order]
            for index in active:
                # XCMS CSV uses ~15 significant digits; include only floating
                # roundtrip error at the original measured mass boundaries.
                epsilon = max(abs(mass_max[index]) * 2e-14, 1e-12)
                left = int(np.searchsorted(masses, mass_min[index] - epsilon, side="left"))
                right = int(np.searchsorted(masses, mass_max[index] + epsilon, side="right"))
                observed = intensities[left:right]
                intensity = float(observed.sum())
                if intensity > maximum[index]:
                    maximum[index] = intensity
                    result.loc[index, ["apex_rt_min", "apex_mz", "apex_intensity", "apex_native_scan_id"]] = (
                        rt, float(np.dot(masses[left:right], observed) / intensity), intensity, str(raw.get("id", "")))
    return result


def match_local_feature(features: pd.DataFrame, spectrum, ppm: float) -> int | None:
    """Require a measured local peak and a qualified measured preceding MS1."""
    refinement = spectrum.metadata.get("precursor_refinement", {})
    if not refinement.get("ms1_mz") or not refinement.get("ms1_native_id"):
        return None
    if refinement.get("status") != "confirmed":
        return None
    parent_rt = refinement["ms1_rt_seconds"] / 60.0
    mz = float(spectrum.precursor_mz)
    error = (features.feature_mz - mz).abs() / features.feature_mz * 1e6
    eligible = features.loc[(error <= ppm) & (features.rt_min <= spectrum.rt_minutes) &
        (spectrum.rt_minutes <= features.rt_max) & (features.rt_min <= parent_rt) &
        (parent_rt <= features.rt_max) & features.apex_rt_min.notna() &
        features.source_path.eq(str(Path(spectrum.metadata["source"]).resolve())) &
        features.polarity.eq(spectrum.polarity)]
    if eligible.empty:
        return None
    # Two overlapping XCMS boundaries do not establish which local peak was
    # fragmented. Avoid assigning one scan twice or inventing a peak identity.
    if len(eligible) != 1:
        return None
    return int(eligible.index[0])


def select_feature_annotations(matches: pd.DataFrame) -> pd.DataFrame:
    """Retain one best acquired spectrum per local peak; keep its native ranking."""
    if matches.empty:
        return matches.copy()
    keys = ["source_path", "polarity", "xcms_peak_id"]
    # Select the acquisition by its strongest native-ranked result. Do not
    # combine fragments, average scores, rerank candidates or mix acquisitions.
    first = matches.loc[matches.groupby(keys + ["native_scan_id"], sort=False, dropna=False).cumcount().eq(0)]
    ordered = first.sort_values(["total_score", "matched_fragment_count"], ascending=[False, False], kind="stable")
    chosen = ordered.drop_duplicates(keys)[keys + ["native_scan_id"]]
    selected = matches.merge(chosen.assign(_chosen=True), on=keys + ["native_scan_id"], how="left", sort=False)
    return selected.loc[selected._chosen.eq(True)].drop(columns="_chosen").reset_index(drop=True)


def run_untargeted(config: UntargetedConfig, logger=None) -> Path:
    """Extract measured MS1 features, identify their acquired MS2, export XLSX."""
    logger = logger or logging.getLogger(__name__)
    for path in config.mzml_files:
        if not path.is_file():
            raise FileNotFoundError(f"mzML 文件不存在：{path}")
    if config.output_file.exists() and not config.output_file.is_file():
        raise ValueError("输出位置须为 XLSX 文件。")
    protected = [library_path("+"), library_path("-")]
    if config.output_file in protected:
        raise ValueError("输出不能覆盖独立鉴定库。")
    find_rscript(config.rscript_path)
    config.output_file.parent.mkdir(parents=True, exist_ok=True)
    frames, all_rows, input_metadata, audits = [], [], [], []
    logger.info("Validating inputs: %s mzML files, no target table", len(config.mzml_files))
    with TemporaryDirectory(prefix=".sphingolipid_xcms_", dir=config.output_file.parent) as temporary:
        with SphingolipidIdentifier(config.identification_config()) as identifier:
            for number, path in enumerate(config.mzml_files, 1):
                work = Path(temporary) / str(number)
                features, versions = run_xcms(path, config, work, logger)
                logger.info("Measuring MS1 peak apices: %s", path.name)
                features = measure_local_apices(features, path, config.polarity)
                features["qualified_ms2_count"] = 0
                features["identified_ms2_count"] = 0
                logger.info("Matching acquired MS2: %s", path.name)
                if not features.empty:
                    for spectrum in iter_prepared_mzml(path):
                        if spectrum.polarity != config.polarity:
                            continue
                        index = match_local_feature(features, spectrum, config.precursor_ppm)
                        refinement = spectrum.metadata.get("precursor_refinement", {})
                        audit = dict(source_file=path.name, source_path=str(path),
                            native_scan_id=spectrum.metadata.get("native_scan_id", ""),
                            ms2_rt_min=spectrum.rt_minutes, precursor_mz=spectrum.precursor_mz,
                            precursor_mz_raw=refinement.get("raw_mz"),
                            preceding_ms1_native_id=refinement.get("ms1_native_id"),
                            ms1_status=refinement.get("status"), xcms_peak_id=None,
                            association_status="no_unique_local_peak" if refinement.get("status") == "confirmed" else "no_confirmed_ms1",
                            passing_candidate_count=0)
                        audits.append(audit)
                        if index is None:
                            continue
                        audit.update(xcms_peak_id=features.loc[index, "xcms_peak_id"], association_status="linked")
                        features.loc[index, "qualified_ms2_count"] += 1
                        rows = identifier.score_spectrum(spectrum)
                        audit["passing_candidate_count"] = len(rows)
                        if rows:
                            features.loc[index, "identified_ms2_count"] += 1
                        for row in rows:
                            row.update(features.loc[index].drop(labels=["qualified_ms2_count", "identified_ms2_count"]).to_dict())
                            row.update(native_scan_id=spectrum.metadata.get("native_scan_id", ""),
                                collision_energy=spectrum.metadata.get("collision_energy"))
                        all_rows.extend(rows)
                features["annotation_status"] = np.where(features.identified_ms2_count > 0, "identified",
                    np.where(features.qualified_ms2_count > 0, "no_passing_ms2_candidate", "no_qualified_local_ms2"))
                frames.append(features)
                entry = {"source_file": path.name, "source_path": str(path), "size_bytes": path.stat().st_size, **versions}
                input_metadata.append(entry)
                logger.info("File completed: %s, %s features, %s identified features", path.name, len(features),
                            int(features.identified_ms2_count.gt(0).sum()))
                if config.save_intermediate:
                    intermediate = config.output_file.parent / (config.output_file.stem + "_xcms") / str(number)
                    intermediate.mkdir(parents=True, exist_ok=True)
                    for artifact in work.iterdir():
                        shutil.copy2(artifact, intermediate / artifact.name)
                    features.to_csv(intermediate / "measured_features.csv", index=False, encoding="utf-8-sig")
    matches = compatible_frame(all_rows)
    if not matches.empty:
        matches["signal_intensity"] = matches["apex_intensity"]
        matches["Abund"] = matches["apex_intensity"]
    annotations = select_feature_annotations(matches)
    features = pd.concat(frames, ignore_index=True)
    settings = {key: str(value) if isinstance(value, Path) else json.dumps(value, default=str) if isinstance(value, tuple) else value
                for key, value in asdict(config).items()}
    from . import __version__
    settings.update(application="SphinGOlipID", version=__version__, feature_detection="XCMS centWave",
        rt_units="minutes (peakwidth is seconds)", feature_grouping="per file; no alignment or gap filling",
        xcms_prefilter_scans=3, xcms_prefilter_intensity=config.noise, xcms_workers=1,
        xcms_mzdiff=-0.001, xcms_mz_center_fun="wMean", xcms_integrate=1, xcms_fitgauss=False,
        centwave_parameters_json=input_metadata[0].get("centwave_parameters_json"),
        ms2_selection="one best acquired spectrum per local feature; native candidate ranks preserved")
    manifest = json.loads(Path(str(files("sphingolipid_toolkit").joinpath("data", "library_revision.json"))).read_text(encoding="utf-8"))
    settings["library_revision_json"] = json.dumps(manifest, sort_keys=True, separators=(",", ":"))
    logger.info("Exporting results: %s", config.output_file)
    write_identification_workbook({"annotations": annotations, "features": features,
        "spectrum_matches": matches, "spectrum_audit": pd.DataFrame(audits, columns=[
            "source_file", "source_path", "native_scan_id", "ms2_rt_min", "precursor_mz", "precursor_mz_raw",
            "preceding_ms1_native_id", "ms1_status", "xcms_peak_id", "association_status", "passing_candidate_count"]),
        "parameters": pd.DataFrame(settings.items(), columns=["parameter", "value"]),
        "input_files": pd.DataFrame(input_metadata)}, config.output_file)
    logger.info("Untargeted analysis completed: %s features; %s annotated features", len(features),
        int(features.identified_ms2_count.gt(0).sum()))
    return config.output_file
