"""Tk review of an exported workbook, with optional ECN and sample areas."""
from __future__ import annotations

from pathlib import Path
from queue import Empty, Queue
from threading import Thread
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import numpy as np
import pandas as pd

from .annotation_report import (DEFAULT_ECN_WINDOW_MIN, assess_ecn, export_annotation_report,
                                initial_review, load_annotation_source)
from .gui_theme import (ACCENT, ALTERNATE, BACKGROUND, INK, MUTED, PANEL, SELECTION,
                        WARNING, style_review_axis)
from .rt_iup import predict_values


PAGE_SIZE = 200
STATUS_NAMES = {'not_run': 'Not reviewed', 'unmodeled': 'Unmodeled',
                'within_window': 'Within window', 'outside_window': 'Outside window'}
TABLE_FIELDS = (
    ('sample_id', 'Sample', 74), ('lipid_name', 'Lipid name', 156), ('apex_rt_min', 'Apex RT', 90),
    ('predicted_rt_min', 'Fit RT', 77), ('delta_rt_min', 'ΔRT', 62), ('ecn_status', 'ECN status', 132),
)


class AnnotationReviewDialog:
    def __init__(self, parent: tk.Misc, workbook: Path):
        self.window = tk.Toplevel(parent)
        self.window.title('SphinGOlipID · ECN Review & Export')
        self.window.configure(background=BACKGROUND)
        self.window.transient(parent)
        self.window.geometry(f'{min(1120, self.window.winfo_screenwidth()-60)}x{min(680, self.window.winfo_screenheight()-100)}')
        self.window.minsize(820, 480)
        self.window.protocol('WM_DELETE_WINDOW', self._close)
        self.window.bind('<Destroy>', self._cancel_poll)
        self.window.grab_set()
        self.review = None
        self.timer = None
        self.busy = False
        self.saving = False
        self.page = 0
        self.visible = pd.DataFrame()
        self.rt = tk.StringVar(self.window, value=f'{DEFAULT_ECN_WINDOW_MIN:.2f}')
        self.exclude = tk.BooleanVar(self.window, value=False)
        self.status = tk.StringVar(self.window)
        self.page_label = tk.StringVar(self.window)
        self._build(workbook)
        self.rt.trace_add('write', self._clear_review)
        self._submit(lambda: initial_review(load_annotation_source(workbook)), self._show_review, 'Reading annotations…')

    def _build(self, workbook: Path):
        body = ttk.Frame(self.window, padding=10)
        body.pack(fill='both', expand=True)
        body.columnconfigure(0, weight=1)
        body.rowconfigure(4, weight=1)
        ttk.Label(body, text='ECN Review & Export', style='Section.TLabel').grid(row=0, column=0, sticky='w')
        ttk.Label(body, text=workbook.name, style='Muted.TLabel').grid(row=1, column=0, sticky='w', pady=(2, 7))

        controls = ttk.Frame(body)
        controls.grid(row=2, column=0, sticky='ew')
        ttk.Label(controls, text='RT window ± (min)').pack(side='left')
        self.rt_entry = ttk.Spinbox(controls, textvariable=self.rt, from_=.01, to=10, increment=.05, width=6)
        self.rt_entry.pack(side='left', padx=(6, 10))
        self.run_button = ttk.Button(controls, text='Run ECN', command=self._run)
        self.run_button.pack(side='left')
        ttk.Label(controls, text='Series').pack(side='left', padx=(18, 6))
        self.series = ttk.Combobox(controls, width=44, values=('All annotations',), state='disabled')
        self.series.current(0)
        self.series.pack(side='left', fill='x', expand=True)
        self.series.bind('<<ComboboxSelected>>', lambda _event: self._select_series())
        ttk.Label(body, textvariable=self.status, style='Muted.TLabel').grid(row=3, column=0, sticky='w', pady=(5, 8))

        panes = ttk.Panedwindow(body, orient='horizontal')
        panes.grid(row=4, column=0, sticky='nsew')
        table = ttk.Frame(panes, width=604)
        panes.add(table, weight=5)
        self.table = ttk.Treeview(table, columns=[key for key, *_ in TABLE_FIELDS], show='headings', selectmode='browse')
        for key, title, width in TABLE_FIELDS:
            self.table.heading(key, text=title)
            self.table.column(key, width=width, anchor='w' if key in {'lipid_name', 'ecn_status'} else 'center', stretch=False)
        self.table.tag_configure('outside_window', foreground=WARNING)
        self.table.tag_configure('alternate', background=ALTERNATE)
        ybar = ttk.Scrollbar(table, orient='vertical', command=self.table.yview)
        xbar = ttk.Scrollbar(table, orient='horizontal', command=self.table.xview)
        self.table.configure(yscrollcommand=ybar.set, xscrollcommand=xbar.set)
        table.columnconfigure(0, weight=1)
        table.rowconfigure(0, weight=1)
        self.table.grid(row=0, column=0, sticky='nsew')
        ybar.grid(row=0, column=1, sticky='ns')
        xbar.grid(row=1, column=0, sticky='ew')
        pages = ttk.Frame(table)
        pages.grid(row=2, column=0, sticky='ew', pady=(5, 0))
        ttk.Label(pages, textvariable=self.page_label, style='Muted.TLabel').pack(side='left')
        self.next_button = ttk.Button(pages, text='Next', command=lambda: self._turn_page(1), state='disabled')
        self.next_button.pack(side='right')
        self.previous_button = ttk.Button(pages, text='Previous', command=lambda: self._turn_page(-1), state='disabled')
        self.previous_button.pack(side='right', padx=(0, 6))

        chart = ttk.Frame(panes, style='Bordered.TFrame', padding=1)
        panes.add(chart, weight=4)
        self.figure = Figure(figsize=(4.25, 3.1), dpi=100, facecolor=PANEL)
        self.canvas = FigureCanvasTkAgg(self.figure, master=chart)
        self.canvas.get_tk_widget().pack(fill='both', expand=True)
        self._draw_plot(None)

        footer = ttk.Frame(body)
        footer.grid(row=5, column=0, sticky='ew', pady=(8, 0))
        self.exclude_box = ttk.Checkbutton(footer, text='Exclude candidates outside the RT window', variable=self.exclude)
        self.exclude_box.pack(side='left')
        self.export_button = ttk.Button(footer, text='Export Annotation Table', command=self._export,
                                       style='Primary.TButton')
        self.export_button.pack(side='right')

    def _refresh_controls(self):
        available = self.review is not None and not self.busy
        fit_possible = available and not self.review.rows.empty
        has_outside = available and self.review.rows.ecn_status.eq('outside_window').any()
        self.run_button.configure(state='normal' if fit_possible else 'disabled')
        self.rt_entry.configure(state='normal' if fit_possible else 'disabled')
        self.series.configure(state='readonly' if available else 'disabled')
        self.exclude_box.configure(state='normal' if has_outside else 'disabled')
        self.export_button.configure(state='normal' if available else 'disabled')

    def _submit(self, work, done, message: str):
        result = Queue()
        self.busy = True
        self.status.set(message)
        self._refresh_controls()

        def execute():
            try:
                result.put((None, work()))
            except Exception as exc:
                result.put((exc, None))

        def poll():
            try:
                error, value = result.get_nowait()
            except Empty:
                self.timer = self.window.after(60, poll)
                return
            self.timer = None
            self.busy = self.saving = False
            if error is not None:
                self._show_summary()
                messagebox.showerror('ECN Review & Export', str(error), parent=self.window)
            else:
                done(value)
            self._refresh_controls()

        self.timer = self.window.after(60, poll)
        Thread(target=execute, daemon=not self.saving).start()

    def _show_summary(self):
        if self.review is None:
            self.status.set('Could not load annotations.')
        elif self.review.rt_window_min is None:
            self.status.set('Run ECN or export the annotation table.' if not self.review.rows.empty else 'No annotations in this workbook.')
        else:
            counts = self.review.rows.ecn_status.value_counts()
            self.status.set('   ·   '.join(f'{counts.get(key, 0)} {STATUS_NAMES[key].lower()}'
                                          for key in ('within_window', 'outside_window', 'unmodeled')))

    def _show_review(self, review):
        self.review = review
        labels = [f"{row.sample_id} · {row.ecn_family.replace('x:y', f'x:{row.acyl_double_bonds}')} · {row.adduct}"
                  for row in review.models.itertuples()]
        self.series.configure(values=('All annotations', *labels))
        self.series.current(1 if labels else 0)
        self.exclude.set(False)
        self._select_series()
        self._show_summary()

    def _clear_review(self, *_args):
        if self.review is not None and self.review.rt_window_min is not None:
            self._show_review(initial_review(self.review.source))
            self._refresh_controls()

    def _run(self):
        if self.busy or self.review is None:
            return
        try:
            window = float(self.rt.get())
        except ValueError:
            messagebox.showerror('RT window', 'Enter the RT window in minutes.', parent=self.window)
            return
        source = self.review.source
        self._submit(lambda: assess_ecn(source, window), self._show_review, 'Fitting ECN series…')

    def _select_series(self):
        model_index = self.series.current() - 1
        self.visible = self.review.rows
        model = None if model_index < 0 else self.review.models.iloc[model_index]
        if model is not None:
            self.visible = self.visible.loc[self.visible.model_id.eq(model.model_id).fillna(False)]
        self.page = 0
        self._show_page()
        self._draw_plot(model)

    def _turn_page(self, direction):
        self.page = max(0, min((len(self.visible) - 1) // PAGE_SIZE, self.page + direction))
        self._show_page()

    def _show_page(self):
        self.table.delete(*self.table.get_children())
        start = self.page * PAGE_SIZE
        items = self.visible.iloc[start:start + PAGE_SIZE]
        for index, row in enumerate(items.itertuples()):
            times = [f'{value:.3f}' if pd.notna(value) else ''
                     for value in (row.apex_rt_min, row.predicted_rt_min, row.delta_rt_min)]
            self.table.insert('', 'end', values=(row.sample_id, row.lipid_name, *times, STATUS_NAMES[row.ecn_status]),
                              tags=(row.ecn_status, 'alternate') if index % 2 else (row.ecn_status,))
        count, last = len(self.visible), start + len(items)
        self.page_label.set(f'{start+1}–{last} of {count}' if count else '0 annotations')
        self.previous_button.configure(state='normal' if self.page else 'disabled')
        self.next_button.configure(state='normal' if last < count else 'disabled')

    def _draw_plot(self, model):
        self.figure.clear()
        axis = self.figure.add_subplot(111)
        if model is None:
            axis.set_axis_off()
            axis.text(.5, .5, 'Select a fitted series', transform=axis.transAxes, ha='center', va='center', color=MUTED)
        else:
            self.figure.subplots_adjust(left=.17, right=.95, top=.94, bottom=.17)
            start, stop = model.fit_carbon_min, model.fit_carbon_max
            x = self.visible.acyl_carbon
            rt = self.visible.apex_rt_min
            interval = np.linspace(min(x.min(), start), max(x.max(), stop), 150)
            window = self.review.rt_window_min
            predicted = predict_values(model.fit_type, model.coefficients, interval)
            axis.fill_between(interval, predicted - window, predicted + window, color=SELECTION, lw=0)
            axis.plot(interval, predicted, color=ACCENT, ls='--', lw=.9)
            fitted = np.linspace(start, stop, 100)
            axis.plot(fitted, predict_values(model.fit_type, model.coefficients, fitted), color=ACCENT, lw=1.35)
            outside = self.visible.ecn_status.eq('outside_window')
            axis.scatter(x[~outside], rt[~outside], color=ACCENT, s=25, zorder=3)
            axis.scatter(x[outside], rt[outside], facecolor=PANEL, edgecolor=WARNING, s=29, zorder=3)
            style_review_axis(axis)
            axis.set_xlabel('Acyl carbon number', fontsize=10, color=INK)
            axis.set_ylabel('RT (min)', fontsize=10, color=INK)
        self.canvas.draw_idle()

    def _export(self):
        if self.busy or self.review is None:
            return
        review, exclude = self.review, self.exclude.get()
        name = filedialog.asksaveasfilename(parent=self.window, title='Export Annotation Table', defaultextension='.xlsx',
                                            filetypes=[('Excel workbook', '*.xlsx')], initialdir=review.source.path.parent,
                                            initialfile=f'{review.source.path.stem}_annotations.xlsx')
        if name:
            self.saving = True
            self._submit(lambda: export_annotation_report(review, Path(name), exclude_outside=exclude),
                         lambda path: self.status.set(f'Exported {path.name}'), 'Writing the Excel workbook…')

    def _close(self):
        if not self.saving:
            self.window.destroy()

    def _cancel_poll(self, event):
        if event.widget is self.window and self.timer is not None:
            self.window.after_cancel(self.timer)
