import os
from pathlib import Path

import pandas as pd
import pytest

from sphingolipid_toolkit.annotation_report import (
    _ecn_coordinates, annotation_report_tables, assess_ecn, export_annotation_report,
    initial_review, load_annotation_source,
)
from sphingolipid_toolkit.result_writer import write_identification_workbook


@pytest.fixture
def result_workbook(tmp_path):
    paths = [str(tmp_path / folder / '=sample.mzML') for folder in ('a', 'b', 'unused')]

    def candidate(file, peak, carbon, mz, rt, area, *, rank=1, name=None):
        return dict(source_path=paths[file], xcms_peak_id=peak, polarity='-', adduct='[M-H]-',
                    lipid_name=name or f'Cer(d18:1/{carbon}:0)', observed_mz=mz, apex_rt_min=rt,
                    peak_area=area, rank=rank, total_score=91.25, resolution_level='chain_level',
                    native_scan_id='00121', ms2_evidence_json='{"evidence":"' + 'x'*33000 + '"}')

    rows = [
        candidate(0, 'CP1', 16, 536.50, 10.0, 105),
        candidate(0, 'CP2', 18, 564.50, 11.0, 210),
        candidate(0, 'CP3', 20, 592.50, 12.0, 317),
        candidate(0, 'CP4', 22, 620.50, 14.4, 448, rank=2),
        candidate(0, 'CP4', 22, 620.50, 14.4, 448, rank=2, name='Cer(d40:1)'),
        candidate(1, 'CP1', 16, 536.5001, 10.12, 0),
        candidate(1, 'CP2', 18, 564.50, 11.6, 221, rank=2),
        candidate(1, 'CP3', 20, 592.50, 12.05, 333, rank=2),
        candidate(1, 'CP4', 20, 592.50, 12.1, 445, rank=2),
    ]
    path = tmp_path / 'original.xlsx'
    return write_identification_workbook({
        'annotations': pd.DataFrame(rows),
        'features': pd.DataFrame({'annotation_status': ['identified']}),
        'input_files': pd.DataFrame({'source_file': ['=sample.mzML']*3, 'source_path': paths}),
        'parameters': pd.DataFrame({'parameter': ['precursor_ppm'], 'value': [10]}),
    }, path)


def test_ecn_applies_local_top1_model_without_reranking_or_inventing_chains(result_workbook):
    source = load_annotation_source(result_workbook)
    review = assess_ecn(source)
    assert review.models.sample_id.tolist() == ['S01']
    assert review.models[['fit_strategy', 'fit_inliers', 'final_inlier_count']].values.tolist() == [['IOLS', 3, 3]]
    assert review.rows.ecn_status.tolist() == (['within_window']*3 + ['outside_window'] + ['unmodeled']*5)
    assert review.rows.fit_input.tolist()[:5] == [True, True, True, False, False]
    assert review.rows.loc[3, 'predicted_rt_min'] == pytest.approx(13)
    assert review.rows.loc[3, 'extrapolated']
    assert review.rows['rank'].tolist() == source.annotations['rank'].tolist()
    assert review.rows.total_score.tolist() == source.annotations.total_score.tolist()
    assert _ecn_coordinates('Cer 18:1;O2/22:0') == _ecn_coordinates('Cer(d18:1/22:0)')
    assert _ecn_coordinates('Cer(d18:1/h22:0)') != _ecn_coordinates('Cer(d18:1/22:0;O)')
    assert _ecn_coordinates('Cer(d40:1)') == (None, None, None)


def test_report_preserves_ambiguous_peaks_zero_areas_and_audit(result_workbook, tmp_path):
    source = load_annotation_source(result_workbook)
    report = annotation_report_tables(initial_review(source))
    first, second, unused = report['samples'].area_column
    area_table = report['identification_table']
    shared = area_table.loc[area_table.lipid_name.eq('Cer(d18:1/16:0)')].iloc[0]
    assert shared[first] == 105 and shared[second] == 0 and pd.isna(shared[unused])
    assert len(area_table.loc[area_table.lipid_name.eq('Cer(d18:1/18:0)')]) == 2
    assert len(area_table.loc[area_table.lipid_name.eq('Cer(d18:1/20:0)')]) == 3
    assert area_table[unused].isna().all()

    reviewed = assess_ecn(source)
    before = result_workbook.read_bytes()
    path = export_annotation_report(reviewed, tmp_path / 'report' / 'table.xlsx', exclude_outside=True)
    assert result_workbook.read_bytes() == before
    before_rows = pd.read_excel(result_workbook, sheet_name='annotations', dtype={'native_scan_id': str})
    copied = pd.read_excel(path, sheet_name='annotations', dtype={'native_scan_id': str})
    pd.testing.assert_frame_equal(copied, before_rows)
    assert copied.native_scan_id.iloc[0] == '00121'
    sidecar = Path(copied.evidence_file.iloc[0])
    assert (path.parent / sidecar).read_bytes() == (result_workbook.parent / sidecar).read_bytes()
    with pd.ExcelFile(path) as workbook:
        assert workbook.sheet_names[0] == 'identification_table'
        assert workbook.book['samples']['B2'].data_type == 's'
        output = pd.read_excel(workbook, sheet_name='identification_table')
        assert 'Cer(d18:1/22:0)' not in output.lipid_name.values
        assert 'Cer(d40:1)' in output.lipid_name.values
        trace = pd.read_excel(workbook, sheet_name='area_sources')
        omitted = trace.loc[~trace.included]
        assert omitted.annotation_excel_row.tolist() == [5] and omitted.report_id.isna().all()


def test_report_rejects_changed_source_and_can_export_an_empty_result(result_workbook, tmp_path):
    source = load_annotation_source(result_workbook)
    stat = result_workbook.stat()
    os.utime(result_workbook, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1_000_000_000))
    target = tmp_path / 'stale.xlsx'
    with pytest.raises(ValueError, match='已更新'):
        export_annotation_report(initial_review(source), target)
    assert not target.exists()

    empty = write_identification_workbook({'annotations': pd.DataFrame(columns=['lipid_name']),
                                            'input_files': source.samples[['source_file', 'source_path']]},
                                           tmp_path / 'empty.xlsx')
    reviewed = assess_ecn(load_annotation_source(empty))
    tables = annotation_report_tables(reviewed)
    assert reviewed.rows.empty and reviewed.models.empty
    assert tables['identification_table'].empty
    assert set(source.samples.area_column) <= set(tables['identification_table'].columns)
