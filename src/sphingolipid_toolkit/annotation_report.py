"""Post-identification ECN assessment and measured sample-area export."""
from __future__ import annotations

from dataclasses import dataclass, replace
import json
import math
from pathlib import Path
import re

import numpy as np
import pandas as pd

from .result_writer import SCAN_ID_DTYPES, copy_identification_bundle, workbook_stamp
from .rt_iup import RTIUPConfig, fit_ecn, predict_values


AREA_RT_SPAN_MIN = 0.20
DEFAULT_ECN_WINDOW_MIN = RTIUPConfig().rt_window_min
ANNOTATION_FIELDS = ('source_path', 'polarity', 'xcms_peak_id', 'lipid_name', 'adduct',
                     'observed_mz', 'apex_rt_min', 'peak_area')
OPTIONAL_FIELDS = ('rank', 'total_score', 'resolution_level', 'native_scan_id')
MODEL_FIELDS = ('model_id', 'sample_id', 'polarity', 'adduct', 'ecn_family', 'acyl_double_bonds',
                'fit_type', 'coefficients', 'r_squared', 'fit_carbon_min', 'fit_carbon_max', 'fit_inliers',
                'fit_strategy', 'final_inlier_count')
SOURCE_FIELDS = ('report_id', 'annotation_excel_row', 'sample_id', 'source_path', 'polarity', 'xcms_peak_id',
                 'native_scan_id', 'lipid_name', 'adduct', 'observed_mz', 'apex_rt_min', 'peak_area',
                 'rank', 'total_score', 'ecn_status', 'included')
_NAME = re.compile(r'([^()]+?)\s*(?:\(([^()]+)\)|\s+([^()]+))')
_BASE = re.compile(r'([mdt]?)(\d{1,3}):(\d{1,2})(;O\d*)?')
_ACYL = re.compile(r'(h?)(\d{1,3}):(\d{1,2})(;O\d*)?')
_BASIS = {'m': 1, 'd': 2, 't': 3}


@dataclass(frozen=True)
class AnnotationSource:
    path: Path
    stamp: tuple[int, int, int]
    annotations: pd.DataFrame
    samples: pd.DataFrame
    precursor_ppm: float


@dataclass(frozen=True)
class ECNReview:
    source: AnnotationSource
    rows: pd.DataFrame
    models: pd.DataFrame
    rt_window_min: float | None = None


def load_annotation_source(path: Path) -> AnnotationSource:
    path = path.resolve()
    stamp = workbook_stamp(path)
    with pd.ExcelFile(path, engine='openpyxl') as workbook:
        if not {'annotations', 'input_files'} <= set(workbook.sheet_names):
            raise ValueError('请选择包含 annotations 和 input_files 的非靶向结果工作簿。')
        options = dict(keep_default_na=False, na_values=[''])
        annotations = workbook.parse('annotations', dtype=SCAN_ID_DTYPES, **options)
        samples = workbook.parse('input_files', usecols=['source_file', 'source_path'], **options)
        settings = workbook.parse('parameters', **options) if 'parameters' in workbook.sheet_names else None
    if workbook_stamp(path) != stamp:
        raise ValueError('来源工作簿正在更新，请重新打开。')
    missing = set(ANNOTATION_FIELDS) - set(annotations)
    if missing and not annotations.empty:
        raise ValueError('工作簿缺少鉴定字段：' + ', '.join(sorted(missing)))
    if samples.empty or samples.isna().any().any() or samples.source_path.duplicated().any():
        raise ValueError('input_files 中的样本文件缺失或重复。')

    # Preserve the workbook's sample and annotation order, including samples
    # with no annotations and measurements with a genuine area of zero.
    samples.insert(0, 'sample_id', [f'S{i:02d}' for i in range(1, len(samples) + 1)])
    samples['area_column'] = samples.sample_id + ' | ' + samples.source_file.astype(str)
    annotations = annotations.reindex(columns=ANNOTATION_FIELDS + OPTIONAL_FIELDS)
    annotations.insert(0, 'annotation_excel_row', np.arange(len(annotations)) + 2)
    annotations['sample_id'] = annotations.source_path.map(samples.set_index('source_path').sample_id)
    keys = ['sample_id', 'xcms_peak_id', 'polarity', 'lipid_name', 'adduct']
    if annotations[keys].isna().any().any():
        raise ValueError('鉴定行缺少样本、色谱峰、极性、脂质名称或加合离子。')
    for column in ('observed_mz', 'apex_rt_min', 'peak_area', 'rank', 'total_score'):
        annotations[column] = pd.to_numeric(annotations[column], errors='coerce')
    areas = annotations.peak_area.dropna()
    if not np.isfinite(areas).all() or areas.lt(0).any():
        raise ValueError('积分峰面积须为非负有限数值。')
    for column in ('observed_mz', 'apex_rt_min', 'rank', 'total_score'):
        annotations[column] = annotations[column].where(np.isfinite(annotations[column]))
    annotations['resolution_level'] = annotations.resolution_level.fillna('')
    tolerance = None if settings is None else settings.set_index('parameter').value.get('precursor_ppm')
    ppm = float(tolerance) if pd.notna(tolerance) else 10.0
    if not math.isfinite(ppm) or ppm <= 0:
        raise ValueError('工作簿的 MS1 质量容差无效。')
    return AnnotationSource(path, stamp, annotations, samples, ppm)


def _ecn_coordinates(name: str) -> tuple[str | None, int | None, int | None]:
    """Use reported explicit chains, with separate acyl oxygen and LCB families."""
    match = _NAME.fullmatch(name.strip())
    if match:
        head, parenthesized, bare = match.groups()
        base, _, acyl = (parenthesized or bare).partition('/')
        lcb, fa = _BASE.fullmatch(base.strip()), _ACYL.fullmatch(acyl.strip())
        if lcb and fa:
            lcb_kind, lc, ldb, oxygen = lcb.groups()
            hydroxyl = int(oxygen[2:] or '1') if oxygen else _BASIS.get(lcb_kind)
            alpha, fc, fdb, extra_oxygen = fa.groups()
            if hydroxyl in (1, 2, 3) and (not lcb_kind or _BASIS[lcb_kind] == hydroxyl) and int(fc) > 0 and int(lc) > 0:
                kind = 'mdt'[hydroxyl - 1]
                suffix = extra_oxygen or ''
                return f'{head.strip()}({kind}{lc}:{ldb}/{alpha}x:y{suffix})', int(fc), int(fdb)
    return None, None, None


def initial_review(source: AnnotationSource) -> ECNReview:
    rows = source.annotations.assign(ecn_status='not_run', ecn_family='', acyl_carbon=np.nan,
                                     acyl_double_bonds=np.nan, model_id=pd.NA, fit_input=False,
                                     predicted_rt_min=np.nan, delta_rt_min=np.nan, extrapolated=False)
    return ECNReview(source, rows, pd.DataFrame(columns=MODEL_FIELDS))


def assess_ecn(source: AnnotationSource, rt_window_min: float = DEFAULT_ECN_WINDOW_MIN) -> ECNReview:
    if not math.isfinite(rt_window_min) or rt_window_min <= 0:
        raise ValueError('ECN 的 RT 窗口须为大于零的有限数值。')
    rows = initial_review(source).rows
    rows['ecn_status'] = 'unmodeled'
    coordinates = [_ecn_coordinates(str(name)) for name in rows.lipid_name]
    terms = ['ecn_family', 'acyl_carbon', 'acyl_double_bonds']
    rows[terms] = pd.DataFrame(coordinates, columns=terms, index=rows.index)
    qualified = rows.loc[rows.ecn_family.notna() & np.isfinite(rows.apex_rt_min) & rows.apex_rt_min.ge(0)]
    config = replace(RTIUPConfig(), rt_window_min=rt_window_min)
    records = []
    # Raw apex times are comparable within their own file and ion channel.
    for (sample, mode, adduct), group in qualified.groupby(['sample_id', 'polarity', 'adduct'], sort=False):
        training = group.loc[group['rank'].eq(1)].drop_duplicates(['xcms_peak_id', *terms])
        if training.empty:
            continue
        rows.loc[training.index, 'fit_input'] = True
        curves = fit_ecn(training.rename(columns={
            'ecn_family': '细类', 'acyl_carbon': 'x碳数', 'acyl_double_bonds': '曲线不饱和度',
            'apex_rt_min': '归一化保留时间',
        }), config, keep_short_series=False).lines
        for curve in curves.to_dict('records'):
            family, double_bonds = curve['细类'], int(curve['不饱和度'])
            hits = group.loc[group.ecn_family.eq(family) & group.acyl_double_bonds.eq(double_bonds)]
            predicted = predict_values(curve['拟合类型'], curve['参数'], hits.acyl_carbon.to_numpy(dtype=float))
            residual = hits.apex_rt_min - predicted
            low, high = curve['x最小'], curve['x最大']
            model_id = len(records) + 1
            rows.loc[hits.index, 'model_id'] = model_id
            rows.loc[hits.index, 'predicted_rt_min'] = predicted
            rows.loc[hits.index, 'delta_rt_min'] = residual
            rows.loc[hits.index, 'extrapolated'] = ~hits.acyl_carbon.between(low, high)
            rows.loc[hits.index, 'ecn_status'] = np.where(residual.abs().le(rt_window_min),
                                                         'within_window', 'outside_window')
            records.append((model_id, sample, mode, adduct, family, double_bonds, curve['拟合类型'],
                            json.dumps(curve['参数']), curve['R²'], low, high, curve['fit_point_count'],
                            curve['fit_strategy'], curve['内点数']))
    return ECNReview(source, rows, pd.DataFrame(records, columns=MODEL_FIELDS), rt_window_min)


def _combine_samples(peaks: pd.DataFrame, precursor_ppm: float) -> bool:
    if len(peaks) <= 1:
        return True
    mz, rt = peaks.observed_mz, peaks.apex_rt_min
    return (not peaks.sample_id.duplicated().any() and np.isfinite(mz).all() and np.isfinite(rt).all()
            and mz.min() > 0 and rt.min() >= 0 and mz.max() - mz.min() <= mz.min() * precursor_ppm * 1e-6
            and rt.max() - rt.min() <= AREA_RT_SPAN_MIN + 1e-12)


def _area_groups(rows: pd.DataFrame, precursor_ppm: float, include: pd.Series) -> list[pd.DataFrame]:
    groups = []
    peak_key = ['sample_id', 'xcms_peak_id']
    identity = ['lipid_name', 'adduct', 'polarity', 'resolution_level']
    for _, lipid in rows.groupby(identity, sort=False):
        peaks = lipid.drop_duplicates(peak_key)
        selections = (lipid,) if _combine_samples(peaks, precursor_ppm) else (part for _, part in lipid.groupby(peak_key, sort=False))
        for group in selections:
            retained = group.loc[include.loc[group.index]]
            if not retained.empty:
                groups.append(retained)
    return sorted(groups, key=lambda group: group.index[0])


def annotation_report_tables(review: ECNReview, *, exclude_outside: bool = False) -> dict[str, pd.DataFrame]:
    if exclude_outside and review.rt_window_min is None:
        raise ValueError('请先运行 ECN 复核。')
    source, rows = review.source, review.rows
    samples = source.samples
    sample_headers = samples.set_index('sample_id').area_column.to_dict()
    include = ~rows.ecn_status.eq('outside_window') if exclude_outside else pd.Series(True, index=rows.index)
    report_ids = pd.Series(pd.NA, dtype='Int64', index=rows.index)
    fields = ('lipid_name', 'adduct', 'polarity', 'resolution_level', 'observed_mz', 'apex_rt_min', 'total_score')
    records = []
    for report_id, group in enumerate(_area_groups(rows, source.precursor_ppm, include), 1):
        report_ids.loc[group.index] = report_id
        # One area per local peak. Original candidate rows and rank remain in
        # the source sheets and the row-to-cell audit.
        if group.groupby('sample_id').peak_area.nunique().gt(1).any():
            raise ValueError('同一样本的同一色谱峰包含不同积分面积。')
        first = group.iloc[0]
        status = group.ecn_status.unique()
        areas = group.groupby('sample_id', sort=False).peak_area.first()
        records.append({'report_id': report_id, **first[list(fields)].to_dict(),
                        'representative_sample': first.sample_id,
                        'ecn_status': status[0] if len(status) == 1 else 'mixed',
                        **{sample_headers[sample]: value for sample, value in areas.items()}})

    exported = pd.DataFrame(records, columns=['report_id', *fields, 'representative_sample', 'ecn_status', *sample_headers.values()])
    provenance = rows.assign(report_id=report_ids, included=include)
    settings = pd.DataFrame({
        'parameter': ['source_workbook', 'area_rt_span_min', 'area_precursor_ppm', 'ecn_rt_window_min', 'exclude_ecn_outside_window'],
        'value': [str(source.path), AREA_RT_SPAN_MIN, source.precursor_ppm, review.rt_window_min, exclude_outside],
    })
    return {'identification_table': exported, 'samples': samples, 'area_sources': provenance[list(SOURCE_FIELDS)],
            'ecn_review': provenance, 'ecn_models': review.models, 'report_settings': settings}


def export_annotation_report(review: ECNReview, path: Path, *, exclude_outside: bool = False) -> Path:
    tables = annotation_report_tables(review, exclude_outside=exclude_outside)
    return copy_identification_bundle(review.source.path, path, report_tables=tables, expected_stamp=review.source.stamp)
