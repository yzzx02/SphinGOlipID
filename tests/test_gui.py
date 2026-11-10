"""Behavior checks for GUI parameter wiring, result review, and run lifecycle."""

import gzip
import json
import threading
import time
import tkinter as tk

import pandas as pd
import pytest

from sphingolipid_toolkit import gui
from sphingolipid_toolkit.config import SphinGOlipIDConfig
from sphingolipid_toolkit.result_writer import write_identifications


def parameter_values(**changes):
    values = {
        key: str(SphinGOlipIDConfig.__dataclass_fields__[key].default)
        for key, *_rest in gui.PARAMETERS
    }
    return values | changes


@pytest.mark.parametrize("changes", [
    {"fragment_ppm": "0"}, {"fragment_ppm": "nan"},
    {"precursor_ppm": "-1"}, {"min_total_score": "101"},
    {"min_total_score": "inf"}, {"top_n": "2.5"},
    {"top_n": "0"}, {"top_n": ""},
])
def test_parameters_reject_invalid_user_input(changes):
    with pytest.raises(ValueError):
        gui.parse_matching_parameters(parameter_values(**changes))


def test_bundled_library_gui_is_ready_with_three_inputs(app,tmp_path):
    select_inputs(app,tmp_path)
    app.selection.ms1_db_file=None
    app._refresh_selection_labels()
    assert app.run_button.instate(['!disabled'])
    config=app._build_config()
    assert config.backend=='lipidgate' and config.ms1_library_path is None


def test_result_filter_uses_literal_names_numeric_scores_and_original_rows():
    frame = pd.DataFrame({
        "lipid_name": ["Cer(d18:1/24:0)", "HexCer(d18:1/24:0)", None, "SM(d18:1/16:0)"],
        "ms2_match_score": ["0.8", "bad", None, "0.9"],
        "rank": ["10", "2", None, "1"],
    }, index=[12, 25, 30, 45])
    before = frame.copy(deep=True)
    filtered = gui.filter_result_rows(frame, "cer(d18:", 50)
    assert list(filtered.index) == [12]
    ordered = gui.filter_result_rows(frame, sort_column="rank")
    assert list(ordered.index) == [45, 25, 12, 30]
    pd.testing.assert_frame_equal(frame, before)


def test_legacy_result_columns_and_empty_filter():
    frame = pd.DataFrame({"注释": ["Cer", "SM"], "匹配度分数": [0.35, 0.8]})
    assert list(gui.filter_result_rows(frame, "cer", 35).index) == [0]
    assert gui.filter_result_rows(frame, "absent").empty
    with pytest.raises(ValueError, match="得分列"):
        gui.filter_result_rows(pd.DataFrame({"name": ["Cer"]}), min_score=50)


@pytest.fixture(scope="module")
def tk_root():
    try:
        root = tk.Tk()
    except tk.TclError as exc:
        pytest.skip(f"Tk display unavailable: {exc}")
    root.withdraw()
    yield root
    root.destroy()


@pytest.fixture
def app(monkeypatch, tk_root):
    root = tk.Toplevel(tk_root)
    root.withdraw()
    instance = gui.SphinGOlipIDApp(root)
    monkeypatch.setattr(gui.messagebox, "showerror", lambda *args, **kwargs: pytest.fail(str(args)))
    monkeypatch.setattr(gui.messagebox, "showwarning", lambda *args, **kwargs: pytest.fail(str(args)))
    yield instance
    instance._on_close()


def result_workbook(path):
    with pd.ExcelWriter(path) as writer:
        pd.DataFrame({
            "lipid_name": ["Cer(d18:1/24:0)", "SM(d18:1/16:0)", "Cer(d18:1/18:0)"],
            "observed_mz": [650.6, 703.5, 566.5],
            "observed_rt": [8.5, 10.2, 5.1],
            "ms2_match_score": [0.6, 0.9, 0.8],
            "rank": [1, 1, 2],
            "matched_fragment_count": [3, 5, 4],
            "extra_evidence": ["first", "second", "third"],
        }).to_excel(writer, index=False, sheet_name="annotations")
        pd.DataFrame({"parameter": ["keep this worksheet"]}).to_excel(writer, index=False, sheet_name="metadata")


def test_sorted_filtered_selection_and_export_preserve_original_workbook(app, tmp_path, monkeypatch):
    source = tmp_path / "annotations.xlsx"
    result_workbook(source)
    assert app._load_result_table(source)
    app.search_var.set("cer(")
    app.score_filter_var.set("70")
    app._sort_results("ms2_match_score")
    assert app.result_tree.get_children() == ("2",)
    app.result_tree.selection_set("2")
    app._on_result_select()
    detail = app.detail_text.get("1.0", tk.END)
    assert "Cer(d18:1/18:0)" in detail
    assert "third" in detail
    destination = tmp_path / "saved.xlsx"
    monkeypatch.setattr(gui.filedialog, "asksaveasfilename", lambda **kwargs: str(destination))
    app._save_current_result_as()
    assert destination.read_bytes() == source.read_bytes()
    assert app.result_summary.get() == "显示 1 / 共 3 条"
    app._reset_filters()
    assert len(app.result_tree.get_children()) == 3


def select_inputs(app, tmp_path):
    for key, name in (("txt_file", "input.txt"), ("precursor_file", "targets.xlsx"), ("ms1_db_file", "library.xlsx")):
        path = tmp_path / name
        path.touch()
        setattr(app.selection, key, path)
    app.selection.output_dir = tmp_path
    app._refresh_selection_labels()


def test_setup_controls_fit_at_minimum_window_size(app):
    app.root.attributes("-alpha", 0)
    app.root.geometry("1000x640")
    app.root.deiconify()
    app.root.update()
    page_height = app.setup_page.winfo_height()
    for widget in app.setup_page.winfo_children():
        assert widget.winfo_y() + widget.winfo_height() <= page_height
    app.root.withdraw()


def test_run_freezes_parameters_preserves_navigation_and_recovers_controls(app, tmp_path, monkeypatch):
    select_inputs(app, tmp_path)
    app.parameter_vars["fragment_ppm"].set("12.5")
    app.parameter_vars["top_n"].set("5")
    app.parameter_vars["precursor_ppm"].set("8")
    app.parameter_vars["min_total_score"].set("55")
    app.save_intermediate.set(True)
    started = threading.Event()
    release = threading.Event()
    snapshots = []

    def backend(config, logger):
        snapshots.append(config)
        logger.info("Reading MS2 file: input.txt")
        started.set()
        assert release.wait(5)
        result = config.output_dir / config.final_result_name
        result_workbook(result)
        return [result]

    monkeypatch.setattr(gui, "run_batch", backend)
    app._start_run()
    assert started.wait(2)
    try:
        assert app.running
        assert app.run_button.instate(["disabled"])
        assert all(control.instate(["disabled"]) for control in app.input_controls)
        app.show_setup()
        app.show_results()
        # Even programmatic changes after start cannot affect the worker snapshot.
        app.parameter_vars["top_n"].set("9")
        assert snapshots[0].top_n == 5
        assert snapshots[0].fragment_ppm == 12.5
        assert snapshots[0].precursor_ppm == 8
        assert snapshots[0].min_total_score == 55
        assert snapshots[0].backend == 'lipidgate'
        assert snapshots[0].save_intermediate
        app.root.update()
    finally:
        release.set()
    deadline = time.monotonic() + 5
    while app.running and time.monotonic() < deadline:
        app.root.update()
        time.sleep(0.02)
    assert not app.running
    assert "解析完成" in app.run_status.get()
    assert len(app.result_tree.get_children()) == 3
    assert "Reading MS2" in app.log_text.get("1.0", tk.END)
    assert app.run_button.instate(["!disabled"])
    assert app.export_button.instate(["!disabled"])


def test_run_failure_restores_controls_and_shows_error(app, tmp_path, monkeypatch):
    select_inputs(app, tmp_path)
    errors = []
    monkeypatch.setattr(gui.messagebox, "showerror", lambda title, text, **kwargs: errors.append(text))
    app.running = True
    app._set_running_ui()
    app.ui_queue.put(("error", "invalid target table"))
    app.root.after_cancel(app.queue_after)
    app._drain_queue()
    assert not app.running
    assert errors == ["invalid target table"]
    assert app.run_button.instate(["!disabled"])
    assert app.empty_title.cget("text") == "解析失败"


def test_filename_braces_reach_backend_as_literal_names(app, tmp_path):
    select_inputs(app, tmp_path)
    txt = tmp_path / "input{sample}.txt"
    txt.touch()
    precursor = tmp_path / "targets{sample}.xlsx"
    precursor.touch()
    app.selection.txt_file = txt
    app.selection.precursor_file = precursor
    config = app._build_config()
    assert config.text_pattern.format(i=1) == txt.name
    assert config.target_pattern.format(i=1) == precursor.name


def test_failed_import_retains_previous_result(app, tmp_path, monkeypatch):
    source = tmp_path / "valid.xlsx"
    result_workbook(source)
    assert app._load_result_table(source)
    errors = []
    monkeypatch.setattr(gui.messagebox, "showerror", lambda *args, **kwargs: errors.append(args))
    assert not app._load_result_table(tmp_path / "missing.xlsx")
    assert errors and app.current_result == source
    assert len(app.result_tree.get_children()) == 3


def test_legacy_result_details_include_fragment_counts(app, tmp_path):
    source = tmp_path / "legacy.xlsx"
    pd.DataFrame({"注释": ["Cer(d18:1/16:0)"], "实际mz": [[100, 200]],
                  "强度": [[50, 60]], "碎片": [["first", "second"]]}).to_excel(source, index=False)
    assert app._load_result_table(source)
    assert app.result_df.matched_fragment_count.tolist() == [2]
    app.result_tree.selection_set("0")
    app._on_result_select()
    details = app.detail_text.get("1.0", tk.END)
    assert "实际mz: [100, 200]" in details
    assert "first" in details and "second" in details


def test_queue_drain_is_bounded_and_remains_scheduled(app):
    for number in range(250):
        app.ui_queue.put(("log", (str(number), "Reading MS2 file")))
    app.root.after_cancel(app.queue_after)
    app._drain_queue()
    assert app.ui_queue.qsize() == 150
    assert app.queue_after in app.root.tk.splitlist(app.root.tk.call("after", "info"))


def test_real_gui_run_processes_spectrum_after_upstream_merge(app, tmp_path, monkeypatch):
    from test_bug_regressions import _write_inputs
    config = _write_inputs(tmp_path)
    app.selection = gui.FileSelection(tmp_path / "hilic-msms-1.txt", tmp_path / "hilic-msms-1.xlsx",
                                      config.ms1_library_path, config.output_dir)
    app._refresh_selection_labels()
    app._start_run()
    deadline = time.monotonic() + 10
    while app.running and time.monotonic() < deadline:
        app.root.update()
        time.sleep(.01)
    assert not app.running
    assert app.current_result == config.output_dir / config.final_result_name
    assert not app.result_df.empty
    assert app.result_df.matched_fragment_count.max() >= 2
    assert app.export_button.instate(["!disabled"])


def test_closing_during_analysis_keeps_export_running(app, monkeypatch):
    warnings = []
    monkeypatch.setattr(gui.messagebox, "showwarning", lambda *args, **kwargs: warnings.append(args))
    app.running = True
    try:
        app._on_close()
        assert warnings and not app.closed
        assert app.root.winfo_exists()
    finally:
        app.running = False


def test_worker_start_failure_restores_gui_controls(app, tmp_path, monkeypatch):
    select_inputs(app, tmp_path)
    errors = []
    def fail_start(self):
        raise RuntimeError("thread unavailable")
    monkeypatch.setattr(gui.threading.Thread, "start", fail_start)
    monkeypatch.setattr(gui.messagebox, "showerror", lambda title, text, **kwargs: errors.append(text))
    app._start_run()
    app.root.after_cancel(app.queue_after)
    app._drain_queue()
    assert not app.running
    assert app.run_button.instate(["!disabled"])
    assert all(control.instate(["!disabled"]) for control in app.input_controls)
    assert "thread unavailable" in errors[0]


def test_current_total_score_uses_same_percent_scale_for_display_and_filter(app, tmp_path):
    source = tmp_path / "current.xlsx"
    frame = pd.DataFrame({"lipid_name": ["Cer A", "Cer B"],
                          "total_score": [43.5, 82.25], "rank": [1, 2],
                          "matched_fragment_count": [3, 5]})
    frame.to_excel(source, index=False)
    original = source.read_bytes()
    assert app._load_result_table(source)
    assert app.score_filter_entry.instate(["!disabled"])
    assert float(app.score_filter_entry.cget("to")) == 100
    app.score_filter_var.set("70")
    app._apply_filters()
    assert app.result_tree.get_children() == ("1",)
    assert app.result_tree.set("1", "total_score") == "82.25"
    assert app.result_df["rank"].tolist() == [1, 2]
    assert source.read_bytes() == original


def test_legacy_scores_display_as_percent_without_rewriting_original(app, tmp_path):
    source = tmp_path / "legacy-score.xlsx"
    result_workbook(source)
    original = source.read_bytes()
    assert app._load_result_table(source)
    assert app.result_tree.set("0", "total_score") == "60.0"
    app.score_filter_var.set("70")
    app._apply_filters()
    assert app.result_tree.get_children() == ("1", "2")
    assert source.read_bytes() == original


@pytest.mark.parametrize("value", ["101", "nan", "not a score"])
def test_invalid_score_filter_keeps_previous_view_and_explains_error(app, tmp_path, value):
    source = tmp_path / "result.xlsx"
    result_workbook(source)
    app._load_result_table(source)
    previous = app.result_tree.get_children()
    app.score_filter_var.set(value)
    app._apply_filters()
    assert app.result_tree.get_children() == previous
    assert "有效" in app.filter_message.get() or "100" in app.filter_message.get()


def test_current_total_score_takes_precedence_over_legacy_alias():
    frame = pd.DataFrame({"total_score": [82.5, 30.0], "ms2_match_score": [0.2, 0.9]})
    original = frame.copy(deep=True)
    assert gui.filter_result_rows(frame, min_score=70).index.tolist() == [0]
    pd.testing.assert_frame_equal(frame, original)


def external_result(tmp_path):
    source_dir = tmp_path / "original"
    source_dir.mkdir()
    source = source_dir / "result.xlsx"
    evidence = json.dumps({"matches": [[264.2686, 99.5]] * 3000}, separators=(",", ":"))
    write_identifications(pd.DataFrame({"lipid_name": ["Cer(d18:1/24:0)"],
        "total_score": [88.1234], "rank": [1], "ms2_evidence_json": [evidence]}), source)
    with pd.ExcelWriter(source, mode="a", engine="openpyxl") as writer:
        pd.DataFrame({"setting": ["original metadata"]}).to_excel(writer, sheet_name="metadata", index=False)
    exported = pd.read_excel(source)
    return source, exported.loc[0, "evidence_file"], exported.loc[0, "evidence_sha256"], evidence


def test_save_as_copies_complete_external_evidence_and_all_sheets(app, tmp_path, monkeypatch):
    source, sidecar, key, evidence = external_result(tmp_path)
    source_bytes = source.read_bytes()
    destination = tmp_path / "saved" / "result.xlsx"
    destination.parent.mkdir()
    assert app._load_result_table(source)
    monkeypatch.setattr(gui.filedialog, "asksaveasfilename", lambda **kwargs: str(destination))
    app._save_current_result_as()
    assert destination.read_bytes() == source_bytes
    result = pd.read_excel(destination)
    evidence_path = destination.parent / result.loc[0, "evidence_file"]
    assert evidence_path.is_file()
    assert json.loads(gzip.decompress(evidence_path.read_bytes()))[key] == evidence
    assert evidence_path.read_bytes() == (source.parent / sidecar).read_bytes()
    assert pd.read_excel(destination, sheet_name="metadata").setting.tolist() == ["original metadata"]
    assert result.total_score.tolist() == [88.1234] and result["rank"].tolist() == [1]


@pytest.mark.parametrize("failure", ["missing", "conflict"])
def test_incomplete_save_as_keeps_existing_workbook_and_evidence(app, tmp_path, monkeypatch, failure):
    source, sidecar, _key, _evidence = external_result(tmp_path)
    destination = tmp_path / "saved" / "result.xlsx"
    destination.parent.mkdir()
    destination.write_bytes(b"previous result")
    if failure == "missing":
        (source.parent / sidecar).unlink()
    else:
        (destination.parent / sidecar).write_bytes(b"different evidence; do not overwrite")
    assert app._load_result_table(source)
    errors = []
    monkeypatch.setattr(gui.messagebox, "showerror", lambda *args, **kwargs: errors.append(args))
    monkeypatch.setattr(gui.filedialog, "asksaveasfilename", lambda **kwargs: str(destination))
    app._save_current_result_as()
    assert errors
    assert destination.read_bytes() == b"previous result"
    if failure == "conflict":
        assert (destination.parent / sidecar).read_bytes() == b"different evidence; do not overwrite"
