"""Shared peak filtering, polarity and conservative MS1 precursor correction."""
from __future__ import annotations

from dataclasses import asdict

from .config import DEFAULT_SEARCH_CONFIG, SearchConfig
from .models import ExperimentalSpectrum, normalize_peaks
from .precursor_refinement import MS1Survey, MS1PrecursorRefiner


def make_refiner(surveys, config: SearchConfig = DEFAULT_SEARCH_CONFIG):
    return MS1PrecursorRefiner(
        surveys, association_ppm=config.ms1_precursor_association_ppm,
        confirmation_ppm=config.ms1_precursor_confirmation_ppm,
        max_gap_seconds=config.ms1_precursor_max_gap_seconds,
    )


def prepare_spectrum(*, scan_id, raw_mz, rt_seconds, raw_peaks, polarity,
                     charge, parent_native, refiner, source, min_relative_intensity):
    peaks = [p for p in normalize_peaks(raw_peaks) if p.relative_intensity >= min_relative_intensity]
    if not peaks:
        return None
    refinement = refiner.refine(float(raw_mz), float(rt_seconds), parent_native)
    return ExperimentalSpectrum(
        scan_id=scan_id, precursor_mz=refinement.matching_mz, rt_minutes=float(rt_seconds) / 60,
        polarity=polarity, peaks=peaks, precursor_charge=charge or None,
        metadata={"source": str(source), "precursor_refinement": asdict(refinement)},
    )


def iter_openms_spectra(experiment, *, source, min_relative_intensity,
                       start_rt_min=None, config=DEFAULT_SEARCH_CONFIG):
    surveys = []
    for index, spectrum in enumerate(experiment):
        if spectrum.getMSLevel() == 1:
            mz, intensity = spectrum.get_peaks()
            surveys.append(MS1Survey(str(spectrum.getNativeID()), f"scan_{index + 1}",
                                     float(spectrum.getRT()), mz, intensity, int(spectrum.getType()) == 1))
    refiner = make_refiner(surveys, config)
    for index, spectrum in enumerate(experiment):
        if spectrum.getMSLevel() != 2 or not spectrum.getPrecursors():
            continue
        if start_rt_min is not None and spectrum.getRT() / 60 < start_rt_min:
            continue
        precursor = spectrum.getPrecursors()[0]
        polarity_code = int(spectrum.getInstrumentSettings().getPolarity())
        polarity = {1: "+", 2: "-"}.get(polarity_code, "")
        parent = str(precursor.getMetaValue("spectrum_ref")) if precursor.metaValueExists("spectrum_ref") else ""
        prepared = prepare_spectrum(
            scan_id=f"scan_{index + 1}", raw_mz=precursor.getMZ(), rt_seconds=spectrum.getRT(),
            raw_peaks=list(zip(*spectrum.get_peaks())), polarity=polarity,
            charge=int(precursor.getCharge()), parent_native=parent, refiner=refiner,
            source=source, min_relative_intensity=min_relative_intensity,
        )
        if prepared is not None:
            yield prepared


def precursor_result_fields(spectrum, theoretical_mz: float) -> dict:
    r = spectrum.metadata.get("precursor_refinement")
    if not r:
        return {}
    return dict(
        precursor_mz_raw=r["raw_mz"], precursor_mz_source=r["source"],
        precursor_refinement_status=r["status"], precursor_ms1_mz=r["ms1_mz"],
        precursor_ms1_scan_id=r["ms1_scan_id"], precursor_ms1_native_id=r["ms1_native_id"],
        precursor_ms1_rt_raw_min=r["ms1_rt_seconds"] / 60 if r["ms1_rt_seconds"] is not None else None,
        precursor_confirmation_count=r["confirmation_count"],
        raw_precursor_ppm_error=(r["raw_mz"] - theoretical_mz) / theoretical_mz * 1e6,
    )
