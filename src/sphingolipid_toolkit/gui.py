"""Tkinter desktop GUI for the SphinGOlipID MS2 workflow.

The layout follows the earlier lab desktop tool while calling the current
project backend instead of the old generated command modules.
"""

from __future__ import annotations

import logging
import queue
import threading
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

from .config import SphinGOlipIDConfig
from .ms2_pipeline import run_batch


APP_TITLE = "鞘脂质谱数据智能解析软件"


@dataclass
class FileSelection:
    txt_file: Path | None = None
    precursor_file: Path | None = None
    ms1_db_file: Path | None = None
    output_dir: Path | None = None


class QueueLogHandler(logging.Handler):
    """Send backend log lines to Tk's UI queue."""

    def __init__(self, ui_queue: "queue.Queue[tuple[str, Any]]") -> None:
        super().__init__()
        self.ui_queue = ui_queue
        self.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s", "%H:%M:%S"))

    def emit(self, record: logging.LogRecord) -> None:
        self.ui_queue.put(("log", self.format(record)))


class SphinGOlipIDApp:
    def __init__(self, root: tk.Tk) -> None:
        self.root = root
        self.root.title(APP_TITLE)
        self.root.geometry("1200x720")
        self.root.minsize(1000, 640)
        self.root.configure(bg="#ffffff")

        self.selection = FileSelection()
        self.ui_queue: "queue.Queue[tuple[str, Any]]" = queue.Queue()
        self.running = False
        self.current_result: Path | None = None
        self.images: dict[str, tk.PhotoImage] = {}

        self._load_images()
        if "ICO_PC.png" in self.images:
            self.root.iconphoto(True, self.images["ICO_PC.png"])

        self.container = tk.Frame(self.root, bg="#ffffff")
        self.container.pack(fill=tk.BOTH, expand=True)
        self.show_login()
        self.root.after(120, self._drain_queue)

    @property
    def asset_dir(self) -> Path:
        return Path(__file__).with_name("gui_assets")

    def _load_images(self) -> None:
        for name in [
            "BTN_B_blue.png",
            "BTN_B_blue_hover.png",
            "BTN_B_red.png",
            "BTN_B_red_hover.png",
            "BTN_L_blue.png",
            "BTN_L_red.png",
            "BTN_L2_blue.png",
            "BTN_L2_red.png",
            "ICO_PC.png",
            "Image1.png",
            "Image2.png",
        ]:
            path = self.asset_dir / name
            if path.exists():
                try:
                    image = tk.PhotoImage(file=str(path))
                    if name == "Image1.png":
                        image = image.subsample(3, 3)
                    self.images[name] = image
                except tk.TclError:
                    pass

    def _clear(self) -> None:
        for child in self.container.winfo_children():
            child.destroy()

    def _legacy_button(
        self,
        parent: tk.Widget,
        text: str,
        command: Any,
        image_name: str = "BTN_L_blue.png",
        font_size: int = 18,
        width: int = 260,
        height: int = 64,
    ) -> tk.Button:
        image = self.images.get(image_name)
        button = tk.Button(
            parent,
            text=text,
            command=command,
            image=image,
            compound=tk.CENTER,
            bg="#ffffff",
            fg="#000000",
            activebackground="#ffffff",
            activeforeground="#ff0000",
            borderwidth=0,
            highlightthickness=0,
            font=("Microsoft YaHei UI", font_size, "bold"),
            cursor="hand2",
            width=width,
            height=height,
        )
        return button

    def _title(self, parent: tk.Widget, text: str, size: int = 30) -> tk.Label:
        return tk.Label(
            parent,
            text=text,
            bg="#ffffff",
            fg="#000000",
            font=("Microsoft YaHei UI", size, "normal"),
        )

    def show_login(self) -> None:
        self._clear()
        frame = tk.Frame(self.container, bg="#ffffff")
        frame.pack(fill=tk.BOTH, expand=True)

        title = self._title(frame, APP_TITLE, 32)
        title.pack(pady=(68, 20))

        logo = self.images.get("Image1.png")
        if logo:
            logo_label = tk.Label(frame, image=logo, bg="#ffffff")
            logo_label.image = logo
            logo_label.pack(pady=(0, 12))

        self._legacy_button(
            frame,
            "进入系统",
            self.show_setup,
            image_name="BTN_B_blue.png",
            font_size=18,
            width=210,
            height=72,
        ).pack(pady=16)

    def show_setup(self) -> None:
        self._clear()
        frame = tk.Frame(self.container, bg="#ffffff")
        frame.pack(fill=tk.BOTH, expand=True)

        header = tk.Frame(frame, bg="#ffffff")
        header.pack(fill=tk.X, padx=56, pady=(38, 12))
        self._title(header, APP_TITLE, 28).pack(side=tk.LEFT)
        small_logo = self.images.get("Image2.png")
        if small_logo:
            logo_label = tk.Label(header, image=small_logo, bg="#ffffff")
            logo_label.image = small_logo
            logo_label.pack(side=tk.RIGHT)

        form = tk.Frame(frame, bg="#ffffff")
        form.pack(expand=True)

        rows = [
            ("TXT数据导入", self._choose_txt, "txt_file"),
            ("XLSX数据导入", self._choose_precursor, "precursor_file"),
            ("MS1数据库导入", self._choose_ms1_db, "ms1_db_file"),
            ("结果存储位置", self._choose_output_dir, "output_dir"),
        ]
        self.status_labels: dict[str, tk.Label] = {}
        for row_index, (label, command, key) in enumerate(rows):
            button = self._legacy_button(form, label, command, image_name="BTN_L_blue.png")
            button.grid(row=row_index, column=0, padx=(0, 24), pady=13, sticky="e")
            status = tk.Label(
                form,
                text="未选择",
                bg="#ffffff",
                fg="#444444",
                anchor="w",
                font=("Microsoft YaHei UI", 13),
                width=44,
            )
            status.grid(row=row_index, column=1, sticky="w")
            self.status_labels[key] = status

        actions = tk.Frame(frame, bg="#ffffff")
        actions.pack(pady=(18, 70))
        self.run_button = self._legacy_button(
            actions,
            "确认",
            self._start_run,
            image_name="BTN_B_red.png",
            font_size=18,
            width=210,
            height=72,
        )
        self.run_button.pack(side=tk.LEFT, padx=18)
        self._legacy_button(
            actions,
            "跳过",
            self.show_results,
            image_name="BTN_B_blue.png",
            font_size=18,
            width=210,
            height=72,
        ).pack(side=tk.LEFT, padx=18)

        self._refresh_selection_labels()

    def show_results(self) -> None:
        self._clear()
        frame = tk.Frame(self.container, bg="#ffffff")
        frame.pack(fill=tk.BOTH, expand=True)

        left = tk.Frame(frame, bg="#ffffff", width=360)
        left.pack(side=tk.LEFT, fill=tk.Y, padx=(22, 12), pady=24)
        left.pack_propagate(False)

        right = tk.Frame(frame, bg="#ffffff")
        right.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=(0, 22), pady=24)

        self._legacy_button(
            left,
            "数据导入",
            self._choose_result_file,
            image_name="BTN_L2_blue.png",
            font_size=16,
            width=270,
            height=56,
        ).pack(pady=(10, 14))
        self._legacy_button(
            left,
            "存储结果",
            self._save_current_result_as,
            image_name="BTN_L2_red.png",
            font_size=16,
            width=270,
            height=56,
        ).pack(pady=(0, 14))
        self._legacy_button(
            left,
            "返回设置",
            self.show_setup,
            image_name="BTN_L2_blue.png",
            font_size=16,
            width=270,
            height=56,
        ).pack(pady=(0, 22))

        log_frame = tk.LabelFrame(
            left,
            text="运行日志",
            bg="#ffffff",
            fg="#000000",
            font=("Microsoft YaHei UI", 12),
        )
        log_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 14))
        self.log_text = tk.Text(log_frame, height=14, bg="#fbfbfb", relief=tk.FLAT, wrap=tk.WORD)
        self.log_text.pack(fill=tk.BOTH, expand=True, padx=8, pady=8)

        detail_frame = tk.LabelFrame(
            left,
            text="结果信息",
            bg="#ffffff",
            fg="#000000",
            font=("Microsoft YaHei UI", 12),
        )
        detail_frame.pack(fill=tk.BOTH, expand=True)
        self.detail_text = tk.Text(detail_frame, height=10, bg="#fbfbfb", relief=tk.FLAT, wrap=tk.WORD)
        self.detail_text.pack(fill=tk.BOTH, expand=True, padx=8, pady=8)

        title_row = tk.Frame(right, bg="#ffffff")
        title_row.pack(fill=tk.X)
        self._title(title_row, "保留时间规律校正后的解析结果", 24).pack(side=tk.LEFT, pady=(0, 18))
        small_logo = self.images.get("Image2.png")
        if small_logo:
            logo_label = tk.Label(title_row, image=small_logo, bg="#ffffff")
            logo_label.image = small_logo
            logo_label.pack(side=tk.RIGHT)

        table_frame = tk.Frame(right, bg="#ffffff")
        table_frame.pack(fill=tk.BOTH, expand=True)
        self.result_tree = ttk.Treeview(table_frame, show="headings")
        y_scroll = ttk.Scrollbar(table_frame, orient=tk.VERTICAL, command=self.result_tree.yview)
        x_scroll = ttk.Scrollbar(table_frame, orient=tk.HORIZONTAL, command=self.result_tree.xview)
        self.result_tree.configure(yscrollcommand=y_scroll.set, xscrollcommand=x_scroll.set)
        self.result_tree.grid(row=0, column=0, sticky="nsew")
        y_scroll.grid(row=0, column=1, sticky="ns")
        x_scroll.grid(row=1, column=0, sticky="ew")
        table_frame.rowconfigure(0, weight=1)
        table_frame.columnconfigure(0, weight=1)
        self.result_tree.bind("<<TreeviewSelect>>", self._on_result_select)

        if self.current_result and self.current_result.exists():
            self._load_result_table(self.current_result)

    def _choose_txt(self) -> None:
        path = filedialog.askopenfilename(title="打开TXT文件", filetypes=[("TXT File", "*.txt"), ("All files", "*")])
        if path:
            self.selection.txt_file = Path(path)
            self._refresh_selection_labels()

    def _choose_precursor(self) -> None:
        path = filedialog.askopenfilename(title="打开XLSX文件", filetypes=[("XLSX File", "*.xlsx"), ("All files", "*")])
        if path:
            self.selection.precursor_file = Path(path)
            self._refresh_selection_labels()

    def _choose_ms1_db(self) -> None:
        path = filedialog.askopenfilename(title="打开MS1数据库文件", filetypes=[("XLSX File", "*.xlsx"), ("All files", "*")])
        if path:
            self.selection.ms1_db_file = Path(path)
            self._refresh_selection_labels()

    def _choose_output_dir(self) -> None:
        path = filedialog.askdirectory(title="设置存储的文件夹")
        if path:
            self.selection.output_dir = Path(path)
            self._refresh_selection_labels()

    def _choose_result_file(self) -> None:
        path = filedialog.askopenfilename(title="打开结果XLSX文件", filetypes=[("XLSX File", "*.xlsx"), ("All files", "*")])
        if path:
            self.current_result = Path(path)
            self._load_result_table(self.current_result)

    def _refresh_selection_labels(self) -> None:
        labels = getattr(self, "status_labels", None)
        if not labels:
            return
        for key, label in labels.items():
            value = getattr(self.selection, key)
            label.configure(text=value.name if isinstance(value, Path) else "未选择")

    def _validate_selection(self) -> str | None:
        missing = []
        if not self.selection.txt_file:
            missing.append("TXT数据")
        if not self.selection.precursor_file:
            missing.append("XLSX数据")
        if not self.selection.ms1_db_file:
            missing.append("MS1数据库")
        if not self.selection.output_dir:
            missing.append("结果存储位置")
        if missing:
            return "请先选择：" + "、".join(missing)
        return None

    def _start_run(self) -> None:
        error = self._validate_selection()
        if error:
            messagebox.showwarning("提示", error)
            return
        if self.running:
            return

        self.running = True
        self.run_button.configure(state=tk.DISABLED, text="运行中")
        self.show_results()
        self._append_log("开始解析，请稍候...")

        worker = threading.Thread(target=self._run_backend, daemon=True)
        worker.start()

    def _run_backend(self) -> None:
        assert self.selection.txt_file is not None
        assert self.selection.precursor_file is not None
        assert self.selection.ms1_db_file is not None
        assert self.selection.output_dir is not None

        logger = logging.getLogger("sphingolipid_toolkit.gui")
        for handler in list(logger.handlers):
            logger.removeHandler(handler)
            handler.close()
        logger.setLevel(logging.INFO)
        logger.propagate = False
        logger.addHandler(QueueLogHandler(self.ui_queue))

        config = SphinGOlipIDConfig(
            raw_ms2_dir=self.selection.txt_file.parent,
            precursor_dir=self.selection.precursor_file.parent,
            ms1_library_path=self.selection.ms1_db_file,
            output_dir=self.selection.output_dir,
            file_indices=(1,),
            text_pattern=self.selection.txt_file.name,
            target_pattern=self.selection.precursor_file.name,
            result_pattern="result.xlsx",
            final_result_name="ms2_annotation_results.xlsx",
            save_intermediate=False,
        )

        try:
            results = run_batch(config, logger=logger)
            final_result = self.selection.output_dir / config.final_result_name
            if not final_result.exists() and results:
                final_result = Path(results[0])
            self.ui_queue.put(("done", final_result))
        except Exception as exc:
            self.ui_queue.put(("error", str(exc)))
        finally:
            for handler in list(logger.handlers):
                logger.removeHandler(handler)
                handler.close()

    def _drain_queue(self) -> None:
        try:
            while True:
                kind, payload = self.ui_queue.get_nowait()
                if kind == "log":
                    self._append_log(str(payload))
                elif kind == "done":
                    self.running = False
                    self.current_result = Path(payload)
                    self._append_log(f"解析完成：{self.current_result}")
                    self._load_result_table(self.current_result)
                    messagebox.showinfo("提示", "解析完成")
                elif kind == "error":
                    self.running = False
                    self._append_log(f"解析失败：{payload}")
                    messagebox.showerror("解析失败", str(payload))
        except queue.Empty:
            pass
        self.root.after(120, self._drain_queue)

    def _append_log(self, text: str) -> None:
        log_widget = getattr(self, "log_text", None)
        if log_widget is None:
            return
        log_widget.insert(tk.END, text + "\n")
        log_widget.see(tk.END)

    def _load_result_table(self, path: Path) -> None:
        if not hasattr(self, "result_tree"):
            return
        try:
            df = pd.read_excel(path)
        except Exception as exc:
            messagebox.showerror("打开失败", str(exc))
            return

        preferred = [
            "lipid_name",
            "target",
            "observed_rt",
            "signal_intensity",
            "matched_fragment_count",
            "ms2_match_score",
            "rank",
            "source_result_file",
        ]
        columns = [col for col in preferred if col in df.columns]
        if not columns:
            columns = [str(col) for col in df.columns[:12]]

        self.result_df = df
        self.result_columns = columns
        self.result_tree.delete(*self.result_tree.get_children())
        self.result_tree.configure(columns=columns)

        for col in columns:
            self.result_tree.heading(col, text=col)
            self.result_tree.column(col, width=140, anchor=tk.CENTER, stretch=True)

        for idx, row in df.iterrows():
            values = [self._format_cell(row.get(col, "")) for col in columns]
            self.result_tree.insert("", tk.END, iid=str(idx), values=values)

        self._set_detail(f"已加载结果文件：{path}\n共 {len(df)} 行，{len(df.columns)} 列。")

    def _on_result_select(self, _event: tk.Event) -> None:
        selected = self.result_tree.selection()
        if not selected or not hasattr(self, "result_df"):
            return
        idx = int(selected[0])
        row = self.result_df.iloc[idx]
        interesting = [
            "lipid_name",
            "target",
            "observed_mz",
            "observed_rt",
            "matched_fragment_count",
            "ms2_match_score",
            "瀹為檯mz",
            "寮哄害",
            "纰庣墖",
        ]
        lines = []
        for col in interesting:
            if col in row.index:
                lines.append(f"{col}: {self._format_cell(row[col])}")
        self._set_detail("\n".join(lines))

    def _set_detail(self, text: str) -> None:
        detail_widget = getattr(self, "detail_text", None)
        if detail_widget is None:
            return
        detail_widget.delete("1.0", tk.END)
        detail_widget.insert(tk.END, text)

    def _save_current_result_as(self) -> None:
        if not self.current_result or not self.current_result.exists():
            messagebox.showwarning("提示", "当前没有可保存的结果文件")
            return
        path = filedialog.asksaveasfilename(
            title="保存XLSX文件",
            filetypes=[("XLSX File", "*.xlsx"), ("All files", "*")],
            defaultextension=".xlsx",
            initialfile=self.current_result.name,
        )
        if not path:
            return
        try:
            df = pd.read_excel(self.current_result)
            df.to_excel(path, index=False)
            messagebox.showinfo("提示", "保存成功")
        except Exception as exc:
            messagebox.showerror("保存失败", str(exc))

    @staticmethod
    def _format_cell(value: Any) -> str:
        if isinstance(value, (list, tuple, dict)):
            text = str(value)
            return text if len(text) <= 120 else text[:117] + "..."
        try:
            if pd.isna(value):
                return ""
        except (TypeError, ValueError):
            pass
        if value is None:
            return ""
        text = str(value)
        return text if len(text) <= 120 else text[:117] + "..."


def main() -> int:
    root = tk.Tk()
    SphinGOlipIDApp(root)
    root.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
