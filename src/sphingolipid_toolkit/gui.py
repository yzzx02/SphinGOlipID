"""Desktop controls and result review for the SphinGOlipID MS2 workflow."""

from __future__ import annotations

import logging
import math
import queue
import shutil
import threading
import tkinter as tk
from dataclasses import dataclass
from pathlib import Path
from time import monotonic
from tkinter import filedialog, messagebox, ttk
from typing import Any, Mapping

import pandas as pd

from .config import SphinGOlipIDConfig
from .ms2_pipeline import run_batch
from .scoring import add_standard_score_columns


APP_TITLE = "鞘脂质谱数据智能解析软件"
BACKGROUND = "#f3f6fb"
INK = "#182b46"
MUTED = "#64748b"
ACCENT = "#2563eb"

# Field, label, spinner bounds, increment, numeric type.
PARAMETERS = (
    ("fragment_ppm", "碎片容差 (ppm)", 0.1, 1000, 1, float),
    ("precursor_ppm", "前体容差 (ppm)", 0.1, 1000, 1, float),
    ("min_total_score", "最低总分 (0–100)", 0, 100, 1, float),
    ("top_n", "每个前体保留 Top-N", 1, 1000, 1, int),
)
NAME_COLUMNS = ("lipid_name", "注释", "name", "annotation", "Lipid")
SCORE_COLUMNS = ("ms2_match_score", "匹配度分数")
RESULT_FIELDS = (
    (("lipid_name", "注释", "name", "annotation"), "脂质名称", 240),
    (("observed_mz", "target"), "前体 m/z", 110),
    (("observed_rt", "RT"), "保留时间 (min)", 130),
    (("total_score", "ms2_match_score", "匹配度分数"), "MS2 得分", 100),
    (("rank",), "候选排名", 90),
    (("matched_fragment_count",), "匹配碎片数", 110),
    (("signal_intensity", "Abund"), "信号强度", 120),
    (("source_result_file",), "来源文件", 180),
)
STAGES = (
    ("Loading MS1", "加载 MS1 数据库"),
    ("Validating inputs", "校验输入文件"),
    ("Reading MS2", "读取原始 MS/MS 数据"),
    ("Extracting precursor", "提取前体离子谱图"),
    ("Searching MS1", "检索 MS1 候选"),
    ("Matching precursor candidates", "匹配前体离子候选"),
    ("Generating theoretical", "生成理论碎片"),
    ("Matching precursor fragment", "匹配 MS2 碎片"),
    ("Exporting results", "写入注释结果"),
)


@dataclass
class FileSelection:
    txt_file: Path | None = None
    precursor_file: Path | None = None
    ms1_db_file: Path | None = None
    output_dir: Path | None = None


def parse_matching_parameters(values: Mapping[str, str]) -> dict[str, float | int]:
    """Validate user-entered values before creating a backend configuration."""
    parsed: dict[str, float | int] = {}
    legacy=(('min_fragment_intensity','最低碎片强度',0,1000000,10,float),
            ('min_matched_fragments','最少匹配碎片数',1,1000,1,int),
            ('min_match_score','最低匹配得分',0,1,.05,float))
    schema = ((PARAMETERS[0],*legacy,PARAMETERS[-1])
              if any(key in values for key,*_ in legacy) else PARAMETERS)
    for key, label, _low, _high, _step, numeric_type in schema:
        try:
            value = numeric_type(values[key].strip())
        except (KeyError, TypeError, ValueError, OverflowError) as exc:
            unit = "整数" if numeric_type is int else "数值"
            raise ValueError(f"{label}须填写有效的{unit}。") from exc
        if not math.isfinite(value):
            raise ValueError(f"{label}须填写有限数值。")
        if key in {"fragment_ppm", "precursor_ppm"} and value <= 0:
            raise ValueError("碎片容差须大于 0 ppm。")
        if key == "min_fragment_intensity" and value < 0:
            raise ValueError("最低碎片强度不能为负数。")
        if key in {"min_matched_fragments", "top_n"} and value < 1:
            raise ValueError(f"{label}须为不小于 1 的整数。")
        if key == "min_match_score" and not 0 <= value <= 1:
            raise ValueError("最低匹配得分须在 0 到 1 之间。")
        if key == 'min_total_score' and not 0 <= value <= 100:
            raise ValueError('最低总分须在 0 到 100 之间。')
        parsed[key] = value
    return parsed


def filter_result_rows(
    frame: pd.DataFrame,
    query: str = "",
    min_score: float | None = None,
    sort_column: str | None = None,
    descending: bool = False,
) -> pd.DataFrame:
    """Filter/sort without changing row identities used by the detail pane."""
    view = frame
    if query.strip():
        columns = [column for column in NAME_COLUMNS if column in view.columns] or list(view.columns)
        matches = pd.Series(False, index=view.index)
        for column in columns:
            matches |= view[column].astype("string").str.contains(
                query.strip(), case=False, regex=False, na=False
            )
        view = view.loc[matches]
    if min_score is not None:
        column = next((col for col in SCORE_COLUMNS if col in view.columns), None)
        if column is None:
            raise ValueError("当前结果没有 MS2 得分列，无法按得分筛选。")
        scores = pd.to_numeric(view[column], errors="coerce")
        view = view.loc[scores >= min_score]
    if sort_column in view.columns:
        values = view[sort_column]
        nonempty = values.notna() & values.astype(str).str.strip().ne("")
        numeric = pd.to_numeric(values, errors="coerce")
        sort_values = numeric if numeric.loc[nonempty].notna().all() else (
            values.astype("string").str.casefold().where(nonempty)
        )
        ordered = sort_values.sort_values(
            ascending=not descending, kind="stable", na_position="last"
        ).index
        view = view.loc[ordered]
    return view


class QueueLogHandler(logging.Handler):
    """Deliver logs and stage information to the Tk main thread."""

    def __init__(self, ui_queue: "queue.Queue[tuple[str, Any]]") -> None:
        super().__init__()
        self.ui_queue = ui_queue
        self.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s", "%H:%M:%S"))

    def emit(self, record: logging.LogRecord) -> None:
        self.ui_queue.put(("log", (self.format(record), record.getMessage())))


class SphinGOlipIDApp:
    def __init__(self, root: tk.Tk) -> None:
        self.root = root
        self.root.title(APP_TITLE)
        self.root.geometry("1180x720")
        self.root.minsize(1000, 640)
        self.root.configure(bg=BACKGROUND)
        self.root.protocol("WM_DELETE_WINDOW", self._on_close)

        self.selection = FileSelection()
        self.ui_queue: "queue.Queue[tuple[str, Any]]" = queue.Queue()
        self.running = False
        self.closed = False
        self.current_result: Path | None = None
        self.result_df = pd.DataFrame()
        self.result_columns: list[str] = []
        self.column_titles: dict[str, str] = {}
        self.column_widths: dict[str, int] = {}
        self.sort_column: str | None = None
        self.sort_descending = False
        self.run_started = 0.0
        self.run_stage = "准备解析"
        self.filter_after: str | None = None
        self.images: dict[str, tk.PhotoImage] = {}
        self.input_controls: list[ttk.Widget] = []

        self.path_vars = {key: tk.StringVar(value="未选择") for key in FileSelection.__dataclass_fields__}
        self.parameter_vars = {
            key: tk.StringVar(value=str(SphinGOlipIDConfig.__dataclass_fields__[key].default))
            for key, *_rest in PARAMETERS
        }
        self.save_intermediate = tk.BooleanVar(value=False)
        self.setup_status = tk.StringVar(value="请选择原始谱图、前体目标表与输出位置；使用内置鞘脂库。")
        self.run_status = tk.StringVar(value="准备就绪 · 选择文件后开始解析，或打开已有结果。")
        self.result_path = tk.StringVar(value="尚未打开结果文件")
        self.result_summary = tk.StringVar(value="尚未加载结果")
        self.search_var = tk.StringVar()
        self.score_filter_var = tk.StringVar()
        self.filter_message = tk.StringVar()

        self._configure_styles()
        self._load_images()
        self._build_shell()
        self._build_setup()
        self._build_results()
        self.search_var.trace_add("write", self._schedule_filter)
        self.score_filter_var.trace_add("write", self._schedule_filter)
        self._refresh_selection_labels()
        self.queue_after = self.root.after(120, self._drain_queue)

    @property
    def asset_dir(self) -> Path:
        return Path(__file__).with_name("gui_assets")

    def _configure_styles(self) -> None:
        style = ttk.Style(self.root)
        style.theme_use("clam")
        style.configure(".", font=("Microsoft YaHei UI", 10))
        style.configure("TFrame", background=BACKGROUND)
        style.configure("TLabel", background=BACKGROUND, foreground=INK)
        style.configure("Muted.TLabel", foreground=MUTED)
        style.configure("Heading.TLabel", font=("Microsoft YaHei UI", 18, "bold"))
        style.configure("Section.TLabel", font=("Microsoft YaHei UI", 13, "bold"))
        style.configure("Card.TLabelframe", background="white", bordercolor="#dce4ef", relief="solid")
        style.configure("Card.TLabelframe.Label", background=BACKGROUND, foreground=INK,
                        font=("Microsoft YaHei UI", 11, "bold"))
        style.configure("Card.TLabel", background="white", foreground=INK)
        style.configure("TButton", padding=(12, 6), background="white", foreground=INK,
                        bordercolor="#d5dfed", focusthickness=1, focuscolor=ACCENT)
        style.map("TButton", background=[("active", "#eaf1ff")], foreground=[("disabled", "#94a3b8")])
        style.configure("Primary.TButton", background=ACCENT, foreground="white", bordercolor=ACCENT)
        style.map("Primary.TButton", background=[("disabled", "#c6d5ee"), ("active", "#1d4ed8")],
                  foreground=[("disabled", "#64748b"), ("!disabled", "white")])
        style.configure("TEntry", padding=5, fieldbackground="white")
        style.configure("TSpinbox", padding=5, fieldbackground="white", arrowsize=14)
        style.configure("TCheckbutton", background="white", foreground=INK)
        style.configure("TNotebook", background=BACKGROUND, borderwidth=0)
        style.configure("TNotebook.Tab", padding=(22, 10), background="#e6edf7", foreground=MUTED)
        style.map("TNotebook.Tab", background=[("selected", "white")], foreground=[("selected", ACCENT)],
                  padding=[("selected", (22, 10)), ("!selected", (22, 10))])
        style.configure("Treeview", background="white", fieldbackground="white", foreground=INK,
                        rowheight=30, bordercolor="#dce4ef")
        style.configure("Treeview.Heading", background="#edf2fa", foreground=INK,
                        padding=(8, 8), font=("Microsoft YaHei UI", 10, "bold"))
        style.map("Treeview", background=[("selected", "#dbeafe")], foreground=[("selected", INK)])
        style.configure("Horizontal.TProgressbar", background=ACCENT, troughcolor="#e2e8f0")

    def _load_images(self) -> None:
        for name in ("ICO_PC.png", "Image2.png"):
            try:
                image = tk.PhotoImage(file=str(self.asset_dir / name))
                if name == "Image2.png":
                    scale = max(1, math.ceil(max(image.width(), image.height()) / 44))
                    image = image.subsample(scale, scale)
                self.images[name] = image
            except tk.TclError:
                continue
        if "ICO_PC.png" in self.images:
            self.root.iconphoto(True, self.images["ICO_PC.png"])

    def _build_shell(self) -> None:
        header = ttk.Frame(self.root, padding=(24, 16, 24, 10))
        header.pack(fill=tk.X)
        if "Image2.png" in self.images:
            ttk.Label(header, image=self.images["Image2.png"]).pack(side=tk.LEFT, padx=(0, 14))
        titles = ttk.Frame(header)
        titles.pack(side=tk.LEFT)
        ttk.Label(titles, text="SphinGOlipID", style="Heading.TLabel").pack(anchor="w")
        ttk.Label(titles, text=APP_TITLE + " · MS2 匹配与注释", style="Muted.TLabel").pack(anchor="w")
        ttk.Label(header, text="数据导入  →  参数设置  →  结果查看", style="Muted.TLabel").pack(side=tk.RIGHT)

        footer = ttk.Frame(self.root, padding=(24, 9))
        footer.pack(side=tk.BOTTOM, fill=tk.X)
        self.progress = ttk.Progressbar(footer, mode="determinate", length=160)
        self.progress.pack(side=tk.RIGHT, padx=(16, 0))
        ttk.Label(footer, textvariable=self.run_status, style="Muted.TLabel").pack(side=tk.LEFT)

        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill=tk.BOTH, expand=True, padx=24, pady=(4, 0))
        self.setup_page = ttk.Frame(self.notebook, padding=(16, 12))
        self.results_page = ttk.Frame(self.notebook, padding=16)
        self.notebook.add(self.setup_page, text="解析设置")
        self.notebook.add(self.results_page, text="MS2 注释结果")

    def _build_setup(self) -> None:
        self.setup_page.columnconfigure(0, weight=1)
        self.setup_page.rowconfigure(4, weight=1)
        inputs = ttk.LabelFrame(self.setup_page, text=" 1 · 输入文件与输出位置 ",
                                style="Card.TLabelframe", padding=(16, 6))
        inputs.grid(row=0, column=0, sticky="ew")
        inputs.columnconfigure(1, weight=1)
        rows = (
            ("原始 MS/MS 数据 (.txt)", "txt_file", self._choose_txt),
            ("前体离子目标表 (.xlsx)", "precursor_file", self._choose_precursor),
            ("结果输出文件夹", "output_dir", self._choose_output_dir),
        )
        for row, (label, key, command) in enumerate(rows):
            ttk.Label(inputs, text=label, style="Card.TLabel").grid(row=row, column=0, sticky="w", padx=(0, 18), pady=3)
            ttk.Entry(inputs, textvariable=self.path_vars[key], state="readonly").grid(row=row, column=1, sticky="ew", pady=3)
            button = ttk.Button(inputs, text="选择文件夹" if key == "output_dir" else "选择文件",
                                command=command, width=10)
            button.grid(row=row, column=2, padx=(12, 0), pady=3)
            self.input_controls.append(button)

        parameters = ttk.LabelFrame(self.setup_page, text=" 2 · MS2 匹配参数 ",
                                    style="Card.TLabelframe", padding=(16, 6))
        parameters.grid(row=1, column=0, sticky="ew", pady=(16, 0))
        for column, (key, label, low, high, step, _numeric_type) in enumerate(PARAMETERS):
            parameters.columnconfigure(column, weight=1, uniform="parameter")
            ttk.Label(parameters, text=label, style="Card.TLabel").grid(row=0, column=column, sticky="w", padx=6, pady=(0, 5))
            spinner = ttk.Spinbox(parameters, textvariable=self.parameter_vars[key],
                                  from_=low, to=high, increment=step, width=10)
            spinner.grid(row=1, column=column, sticky="ew", padx=6)
            self.input_controls.append(spinner)
        check = ttk.Checkbutton(parameters, text="保留中间文件", variable=self.save_intermediate)
        check.grid(row=2, column=0, columnspan=2, sticky="w", padx=6, pady=(6, 0))
        self.input_controls.append(check)
        reset = ttk.Button(parameters, text="恢复默认参数", command=self._reset_parameters)
        reset.grid(row=2, column=3, columnspan=2, sticky="e", padx=6, pady=(6, 0))
        self.input_controls.append(reset)

        actions = ttk.Frame(self.setup_page)
        actions.grid(row=2, column=0, sticky="ew", pady=(10, 0))
        self.run_button = ttk.Button(actions, text="开始解析", command=self._start_run, style="Primary.TButton")
        self.run_button.pack(side=tk.LEFT)
        self.setup_open_button = ttk.Button(actions, text="打开已有结果", command=self._choose_result_file)
        self.setup_open_button.pack(side=tk.LEFT, padx=12)
        ttk.Label(actions, textvariable=self.setup_status, style="Muted.TLabel").pack(side=tk.LEFT, padx=8)

        ttk.Label(self.setup_page, text="输入说明：TXT 为质谱导出的原始谱图；目标表需包含 target、下限、上限、RT、Abund 列。",
                  style="Muted.TLabel", wraplength=850).grid(row=3, column=0, sticky="w", pady=(8, 0))

    def _build_results(self) -> None:
        page = self.results_page
        top = ttk.Frame(page)
        top.pack(fill=tk.X)
        ttk.Label(top, text="MS2 注释结果", style="Section.TLabel").pack(side=tk.LEFT)
        self.export_button = ttk.Button(top, text="另存完整结果", command=self._save_current_result_as, state="disabled")
        self.export_button.pack(side=tk.RIGHT)
        self.results_open_button = ttk.Button(top, text="打开结果文件", command=self._choose_result_file)
        self.results_open_button.pack(side=tk.RIGHT, padx=10)
        ttk.Entry(page, textvariable=self.result_path, state="readonly").pack(fill=tk.X, pady=(12, 10))

        filters = ttk.Frame(page)
        filters.pack(fill=tk.X)
        ttk.Label(filters, text="脂质名称").pack(side=tk.LEFT)
        self.search_entry = ttk.Entry(filters, textvariable=self.search_var, width=25)
        self.search_entry.pack(side=tk.LEFT, padx=(8, 18))
        ttk.Label(filters, text="最低得分").pack(side=tk.LEFT)
        self.score_filter_entry = ttk.Spinbox(filters, textvariable=self.score_filter_var,
                                             from_=0, to=1, increment=0.05, width=7, state="disabled")
        self.score_filter_entry.pack(side=tk.LEFT, padx=8)
        ttk.Button(filters, text="重置筛选", command=self._reset_filters).pack(side=tk.LEFT, padx=8)
        ttk.Label(filters, textvariable=self.result_summary, style="Muted.TLabel").pack(side=tk.RIGHT)
        ttk.Label(page, textvariable=self.filter_message, style="Muted.TLabel").pack(anchor="w", pady=(5, 6))

        panes = ttk.Panedwindow(page, orient=tk.HORIZONTAL)
        panes.pack(fill=tk.BOTH, expand=True)
        table_frame = ttk.Frame(panes)
        sidebar = ttk.Frame(panes, width=270)
        panes.add(table_frame, weight=4)
        panes.add(sidebar, weight=1)
        table_frame.rowconfigure(0, weight=1)
        table_frame.columnconfigure(0, weight=1)
        self.result_tree = ttk.Treeview(table_frame, show="headings", selectmode="browse")
        self.result_tree.tag_configure("alternate", background="#f6f9fe")
        y_scroll = ttk.Scrollbar(table_frame, orient=tk.VERTICAL, command=self.result_tree.yview)
        x_scroll = ttk.Scrollbar(table_frame, orient=tk.HORIZONTAL, command=self.result_tree.xview)
        self.result_tree.configure(yscrollcommand=y_scroll.set, xscrollcommand=x_scroll.set)
        self.result_tree.grid(row=0, column=0, sticky="nsew")
        y_scroll.grid(row=0, column=1, sticky="ns")
        x_scroll.grid(row=1, column=0, sticky="ew")
        self.result_tree.bind("<<TreeviewSelect>>", self._on_result_select)
        self.empty_panel = ttk.Frame(table_frame, padding=22)
        self.empty_title = ttk.Label(self.empty_panel, text="尚未加载注释结果", style="Section.TLabel")
        self.empty_title.pack(pady=(0, 8))
        self.empty_hint = ttk.Label(self.empty_panel, text="完成解析后将在这里显示结果，\n也可以打开已有的 XLSX 结果文件。",
                                    style="Muted.TLabel", justify=tk.CENTER)
        self.empty_hint.pack()
        self.empty_panel.place(relx=0.5, rely=0.45, anchor="center")

        sidebar.columnconfigure(0, weight=1)
        sidebar.rowconfigure(0, weight=3)
        sidebar.rowconfigure(1, weight=2)
        self.detail_text = self._text_panel(sidebar, " 所选候选详情 ", row=0, height=9)
        self.log_text = self._text_panel(sidebar, " 运行日志 ", row=1, height=7)
        self._set_detail("点击表格中的候选行，查看完整注释字段。")

    @staticmethod
    def _text_panel(parent: ttk.Frame, title: str, row: int, height: int) -> tk.Text:
        panel = ttk.LabelFrame(parent, text=title, style="Card.TLabelframe", padding=8)
        panel.grid(row=row, column=0, sticky="nsew", padx=(12, 0), pady=(0 if row == 0 else 12, 0))
        text = tk.Text(panel, width=26, height=height, bg="white", fg=INK, relief=tk.FLAT,
                       wrap=tk.WORD, font=("Microsoft YaHei UI", 9), state=tk.DISABLED)
        scroll = ttk.Scrollbar(panel, orient=tk.VERTICAL, command=text.yview)
        text.configure(yscrollcommand=scroll.set)
        scroll.pack(side=tk.RIGHT, fill=tk.Y)
        text.pack(fill=tk.BOTH, expand=True)
        return text

    def show_login(self) -> None:
        """Compatibility entry point; setup is the initial page."""
        self.show_setup()

    def show_setup(self) -> None:
        self.notebook.select(self.setup_page)

    def show_results(self) -> None:
        self.notebook.select(self.results_page)

    def _choose_file(self, key: str, title: str, extension: str) -> None:
        if self.running:
            return
        path = filedialog.askopenfilename(parent=self.root, title=title,
                                          filetypes=[(extension.upper() + " 文件", "*." + extension)])
        if path:
            setattr(self.selection, key, Path(path))
            self._refresh_selection_labels()

    def _choose_txt(self) -> None:
        self._choose_file("txt_file", "选择原始 MS/MS 数据", "txt")

    def _choose_precursor(self) -> None:
        self._choose_file("precursor_file", "选择前体离子目标表", "xlsx")

    def _choose_ms1_db(self) -> None:
        self._choose_file("ms1_db_file", "选择 MS1 理论数据库", "xlsx")

    def _choose_output_dir(self) -> None:
        if self.running:
            return
        path = filedialog.askdirectory(parent=self.root, title="选择结果输出文件夹")
        if path:
            self.selection.output_dir = Path(path)
            self._refresh_selection_labels()

    def _choose_result_file(self) -> None:
        if self.running:
            return
        path = filedialog.askopenfilename(parent=self.root, title="打开 MS2 注释结果",
                                          filetypes=[("XLSX 结果文件", "*.xlsx")])
        if path and self._load_result_table(Path(path)):
            self.show_results()
            self.run_status.set("已打开结果文件 · 可搜索、筛选或点击表头排序。")

    def _refresh_selection_labels(self) -> None:
        count = 0
        for key, variable in self.path_vars.items():
            path = getattr(self.selection, key)
            variable.set(str(path) if path else "未选择")
            count += path is not None and key in {'txt_file','precursor_file','output_dir'}
        self.setup_status.set(f"已选择 {count}/3 项" if count < 3 else "输入已就绪，使用内置鞘脂库。")
        self.run_button.configure(state=tk.NORMAL if count == 3 and not self.running else tk.DISABLED)

    def _reset_parameters(self) -> None:
        if self.running:
            return
        for key, variable in self.parameter_vars.items():
            variable.set(str(SphinGOlipIDConfig.__dataclass_fields__[key].default))
        self.save_intermediate.set(False)

    def _validate_selection(self) -> str | None:
        labels = {"txt_file": "原始 MS/MS 数据", "precursor_file": "前体离子目标表",
                  "output_dir": "结果输出文件夹"}
        missing = [label for key, label in labels.items() if getattr(self.selection, key) is None]
        if missing:
            return "请先选择：" + "、".join(missing)
        for key in ("txt_file", "precursor_file"):
            if not getattr(self.selection, key).is_file():
                return f"{labels[key]}不存在或不是文件，请重新选择。"
        if self.selection.output_dir.exists() and not self.selection.output_dir.is_dir():
            return "结果输出位置须为文件夹。"
        return None

    def _build_config(self) -> SphinGOlipIDConfig:
        error = self._validate_selection()
        if error:
            raise ValueError(error)
        parameters = parse_matching_parameters({key: var.get() for key, var in self.parameter_vars.items()})
        return SphinGOlipIDConfig(
            raw_ms2_dir=self.selection.txt_file.parent,
            precursor_dir=self.selection.precursor_file.parent,
            ms1_library_path=self.selection.ms1_db_file,
            output_dir=self.selection.output_dir,
            file_indices=(1,),
            text_pattern=self.selection.txt_file.name.replace("{", "{{").replace("}", "}}"),
            target_pattern=self.selection.precursor_file.name.replace("{", "{{").replace("}", "}}"),
            result_pattern="result.xlsx",
            final_result_name="ms2_annotation_results.xlsx",
            save_intermediate=self.save_intermediate.get(),
            enable_rt_validation=False, enable_deduplication=False,
            **parameters,
        )

    def _start_run(self) -> None:
        if self.running:
            return
        try:
            config = self._build_config()
        except ValueError as exc:
            messagebox.showwarning("检查输入与参数", str(exc), parent=self.root)
            return
        self.running = True
        self.run_started = monotonic()
        self.run_stage = "准备解析"
        self.current_result = None
        self._set_running_ui()
        self.score_filter_entry.configure(state=tk.DISABLED)
        self.progress.configure(mode="indeterminate")
        self.progress.start(12)
        self._replace_text(self.log_text, "")
        self.result_df = pd.DataFrame()
        self.result_columns = []
        self.result_path.set("正在生成新的 MS2 注释结果…")
        self._reset_filters()
        self.empty_title.configure(text="正在进行 MS2 匹配")
        self.empty_hint.configure(text="处理完成后自动显示注释结果。\n当前阶段与耗时显示在窗口底部。")
        self._set_detail("正在解析，完成后可选择候选查看详情。")
        self.show_results()
        self._append_log("开始解析。")
        try:
            threading.Thread(target=self._run_backend, args=(config,), daemon=True).start()
        except Exception as exc:
            self.ui_queue.put(("error", f"无法启动后台解析：{exc}"))

    def _run_backend(self, config: SphinGOlipIDConfig) -> None:
        """Use a frozen config snapshot; the worker never reads Tk variables."""
        logger = logging.getLogger(f"sphingolipid_toolkit.gui.{id(self)}")
        logger.setLevel(logging.INFO)
        logger.propagate = False
        handler = QueueLogHandler(self.ui_queue)
        logger.addHandler(handler)
        try:
            results = run_batch(config, logger=logger)
            final_result = config.output_dir / config.final_result_name
            if not final_result.exists() and results:
                final_result = Path(results[0])
            self.ui_queue.put(("done", final_result))
        except Exception as exc:
            self.ui_queue.put(("error", str(exc)))
        finally:
            logger.removeHandler(handler)
            handler.close()

    def _set_running_ui(self) -> None:
        state = tk.DISABLED if self.running else tk.NORMAL
        for widget in self.input_controls:
            widget.configure(state=state)
        self.setup_open_button.configure(state=state)
        self.results_open_button.configure(state=state)
        self.export_button.configure(state=tk.NORMAL if self.current_result and not self.running else tk.DISABLED)
        self.run_button.configure(text="解析中…" if self.running else "开始解析")
        self._refresh_selection_labels()

    def _finish_run(self) -> None:
        self.running = False
        self.progress.stop()
        self.progress.configure(mode="determinate", value=0)
        self._set_running_ui()

    def _drain_queue(self) -> None:
        if self.closed:
            return
        try:
            for _ in range(100):
                kind, payload = self.ui_queue.get_nowait()
                if kind == "log":
                    line, message = payload
                    self._append_log(line)
                    for prefix, stage in STAGES:
                        if message.startswith(prefix):
                            self.run_stage = stage
                            break
                elif kind == "done":
                    self._finish_run()
                    self._append_log(f"解析完成：{payload}")
                    elapsed = int(monotonic() - self.run_started)
                    if self._load_result_table(Path(payload)):
                        self.run_status.set(f"解析完成 · {len(self.result_df)} 条候选 · 耗时 {elapsed} 秒")
                    else:
                        self.run_status.set("解析已结束，但结果文件加载失败。")
                elif kind == "error":
                    self._finish_run()
                    self._append_log(f"解析失败：{payload}")
                    self.run_status.set("解析失败 · 查看运行日志，检查文件与参数后重试。")
                    self.result_path.set("本次解析失败，尚未生成可查看的结果")
                    self.empty_title.configure(text="解析失败")
                    self.empty_hint.configure(text="检查右侧日志，返回解析设置修改后重试。")
                    messagebox.showerror("解析失败", str(payload), parent=self.root)
        except queue.Empty:
            pass
        finally:
            if not self.closed:
                if self.running:
                    seconds = int(monotonic() - self.run_started)
                    self.run_status.set(f"解析中 · {self.run_stage} · 已用时 {seconds // 60:02d}:{seconds % 60:02d}")
                self.queue_after = self.root.after(120, self._drain_queue)

    def _append_log(self, text: str) -> None:
        self.log_text.configure(state=tk.NORMAL)
        self.log_text.insert(tk.END, text + "\n")
        self.log_text.see(tk.END)
        self.log_text.configure(state=tk.DISABLED)

    def _load_result_table(self, path: Path) -> bool:
        try:
            frame = add_standard_score_columns(pd.read_excel(path))
        except Exception as exc:
            messagebox.showerror("打开失败", str(exc), parent=self.root)
            return False
        frame.columns = frame.columns.map(str)
        self.result_df = frame.reset_index(drop=True)
        self.current_result = path
        self.result_path.set(str(path))
        self.result_columns = []
        self.column_titles = {}
        self.column_widths = {}
        for aliases, title, width in RESULT_FIELDS:
            column = next((col for col in aliases if col in frame.columns), None)
            if column is not None:
                self.result_columns.append(column)
                self.column_titles[column] = title
                self.column_widths[column] = width
        if not self.result_columns:
            self.result_columns = list(frame.columns[:12])
            self.column_titles = {col: col for col in self.result_columns}
            self.column_widths = {col: 140 for col in self.result_columns}
        self._reset_filters()
        self.score_filter_entry.configure(state=tk.NORMAL if any(col in frame for col in SCORE_COLUMNS) else tk.DISABLED)
        self.export_button.configure(state=tk.DISABLED if self.running else tk.NORMAL)
        self._set_detail(f"已加载 {len(frame)} 条候选，{len(frame.columns)} 个字段。\n\n点击候选行查看详情。\n点击表头可切换升序／降序。\n\n另存完整结果会保留原文件的全部工作表。")
        return True

    def _schedule_filter(self, *_args: Any) -> None:
        if self.filter_after is not None:
            self.root.after_cancel(self.filter_after)
        self.filter_after = self.root.after(220, self._apply_filters)

    def _reset_filters(self) -> None:
        self.search_var.set("")
        self.score_filter_var.set("")
        self.sort_column = None
        self.sort_descending = False
        self._apply_filters()

    def _sort_results(self, column: str) -> None:
        self.sort_descending = not self.sort_descending if self.sort_column == column else False
        self.sort_column = column
        self._apply_filters()

    def _apply_filters(self) -> None:
        if self.filter_after is not None:
            self.root.after_cancel(self.filter_after)
            self.filter_after = None
        try:
            text = self.score_filter_var.get().strip()
            try:
                minimum = float(text) if text else None
            except ValueError as exc:
                raise ValueError("请输入有效的最低得分，留空显示全部得分。") from exc
            if minimum is not None and (not math.isfinite(minimum) or not 0 <= minimum <= 1):
                raise ValueError("筛选得分须在 0 到 1 之间，留空显示全部得分。")
            view = filter_result_rows(self.result_df, self.search_var.get(), minimum,
                                      self.sort_column, self.sort_descending)
        except ValueError as exc:
            self.filter_message.set(str(exc))
            return
        self.filter_message.set("点击表头排序；最低得分留空时显示全部得分。")
        self.result_tree.delete(*self.result_tree.get_children())
        self.result_tree.configure(columns=self.result_columns)
        for col in self.result_columns:
            suffix = (" ▼" if self.sort_descending else " ▲") if col == self.sort_column else ""
            self.result_tree.heading(col, text=self.column_titles[col] + suffix,
                                     command=lambda column=col: self._sort_results(column))
            self.result_tree.column(col, width=self.column_widths[col], minwidth=75, anchor=tk.W,
                                    stretch=col == self.result_columns[0])
        for display_index, (idx, row) in enumerate(view.iterrows()):
            self.result_tree.insert("", tk.END, iid=str(idx),
                                    values=[self._format_cell(row[col]) for col in self.result_columns],
                                    tags=("alternate",) if display_index % 2 else ())
        self.result_summary.set(f"显示 {len(view)} / 共 {len(self.result_df)} 条")
        if view.empty:
            if self.current_result:
                self.empty_title.configure(text="结果文件中没有候选" if self.result_df.empty else "没有符合筛选条件的候选")
                self.empty_hint.configure(text="尝试调整匹配参数后重新解析。" if self.result_df.empty else "修改脂质名称或最低得分，或重置筛选。")
            self.empty_panel.place(relx=0.5, rely=0.45, anchor="center")
        else:
            self.empty_panel.place_forget()
        self._set_detail("点击表格中的候选行，查看完整注释字段。")

    def _on_result_select(self, _event: tk.Event | None = None) -> None:
        selected = self.result_tree.selection()
        if not selected:
            return
        row = self.result_df.loc[int(selected[0])]
        lines = [f"{self.column_titles.get(str(col), col)}: {self._format_cell(value, limit=2000)}"
                 for col, value in row.items()]
        self._set_detail("\n\n".join(lines))

    @staticmethod
    def _replace_text(widget: tk.Text, text: str) -> None:
        widget.configure(state=tk.NORMAL)
        widget.delete("1.0", tk.END)
        widget.insert(tk.END, text)
        widget.configure(state=tk.DISABLED)

    def _set_detail(self, text: str) -> None:
        self._replace_text(self.detail_text, text)

    def _save_current_result_as(self) -> None:
        if self.running or not self.current_result:
            return
        path = filedialog.asksaveasfilename(parent=self.root, title="另存完整结果（全部工作表）",
                                            filetypes=[("XLSX 结果文件", "*.xlsx")],
                                            defaultextension=".xlsx", initialfile=self.current_result.name)
        if not path:
            return
        try:
            if Path(path).resolve() != self.current_result.resolve():
                shutil.copy2(self.current_result, path)
            self.run_status.set(f"完整结果已另存到：{path}")
        except Exception as exc:
            messagebox.showerror("保存失败", str(exc), parent=self.root)

    @staticmethod
    def _format_cell(value: Any, limit: int = 120) -> str:
        if value is None:
            return ""
        if not isinstance(value, (list, tuple, dict)):
            try:
                if pd.isna(value):
                    return ""
            except (TypeError, ValueError):
                pass
        text = str(value)
        return text if len(text) <= limit else text[:limit - 3] + "..."

    def _on_close(self) -> None:
        if self.closed:
            return
        if self.running:
            messagebox.showwarning("正在解析", "请等待解析完成后再关闭窗口，避免中断结果写入。", parent=self.root)
            return
        self.closed = True
        self.root.after_cancel(self.queue_after)
        if self.filter_after is not None:
            self.root.after_cancel(self.filter_after)
        self.root.destroy()


def main() -> int:
    root = tk.Tk()
    SphinGOlipIDApp(root)
    root.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
