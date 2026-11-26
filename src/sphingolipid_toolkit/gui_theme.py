"""Shared presentation for the Tk desktop windows and their review plots."""

from tkinter import font as tkfont, ttk


BACKGROUND = "#F2F3F5"
PANEL = "#FFFFFF"
HEADER = "#EBEEF1"
BORDER = "#C9CFD6"
RULE = "#E1E5E9"
INK = "#232B33"
MUTED = "#53616D"
DISABLED = "#727D88"
AXIS = "#98A2AC"
ACCENT = "#45698C"
SELECTION = "#DCE9F5"
ALTERNATE = "#F7F8FA"
WARNING = "#A04B22"
UI_FONT = ("Segoe UI", 9)

_BUTTON = "#F1F3F5"
_BUTTON_HOVER = "#E8EDF2"
_BUTTON_PRESSED = "#DCE4EB"
_ACCENT_HOVER = "#3B5C7D"
_ACCENT_PRESSED = "#304D69"
_FOCUS = "#7F9DBA"
_THUMB = "#CFD6DD"


def configure_desktop_style(root):
    style = ttk.Style(root)
    style.theme_use("clam")
    root.configure(background=BACKGROUND)
    style.configure(".", font=UI_FONT, background=BACKGROUND, foreground=INK,
                    bordercolor=BORDER, lightcolor=BORDER, darkcolor=BORDER,
                    focuscolor=ACCENT, troughcolor=HEADER,
                    selectbackground=SELECTION, selectforeground=INK)
    style.map(".", foreground=[("disabled", DISABLED)])
    style.configure("Panel.TFrame", background=PANEL)
    style.configure("Bordered.TFrame", background=BORDER)
    style.configure("TLabel", background=BACKGROUND, foreground=INK)
    style.configure("Muted.TLabel", foreground=MUTED)
    style.configure("Warning.TLabel", foreground=WARNING)
    style.configure("Heading.TLabel", font=(UI_FONT[0], 14, "bold"))
    style.configure("Section.TLabel", font=(UI_FONT[0], 11, "bold"))
    style.configure("Panel.TLabelframe", background=PANEL, borderwidth=1, relief="solid")
    style.configure("Panel.TLabelframe.Label", background=PANEL, foreground=INK,
                    font=(UI_FONT[0], 10, "bold"))
    style.configure("Panel.TLabel", background=PANEL)
    style.configure("TSeparator", background=BORDER)

    style.configure("TButton", padding=(10, 4), background=_BUTTON, foreground=INK,
                    relief="solid", borderwidth=1, focusthickness=1)
    style.map("TButton",
              background=[("disabled", BACKGROUND), ("pressed", _BUTTON_PRESSED), ("active", _BUTTON_HOVER)],
              bordercolor=[("focus", _FOCUS)], foreground=[("disabled", DISABLED)])
    style.configure("Primary.TButton", background=ACCENT, foreground=PANEL,
                    bordercolor=ACCENT, lightcolor=ACCENT, darkcolor=ACCENT, focuscolor=PANEL)
    style.map("Primary.TButton",
              background=[("disabled", HEADER), ("pressed", _ACCENT_PRESSED), ("active", _ACCENT_HOVER)],
              bordercolor=[("disabled", BORDER), ("focus", _ACCENT_PRESSED)],
              lightcolor=[("disabled", BORDER)], darkcolor=[("disabled", BORDER)],
              foreground=[("disabled", DISABLED), ("!disabled", PANEL)])

    for field in ("TEntry", "TSpinbox", "TCombobox"):
        style.configure(field, fieldbackground=PANEL, foreground=INK, insertcolor=INK,
                        borderwidth=1, padding=(5, 3), relief="flat")
        style.map(field, fieldbackground=[("disabled", BACKGROUND), ("readonly", PANEL)],
                  foreground=[("disabled", DISABLED)], bordercolor=[("focus", _FOCUS)],
                  lightcolor=[("focus", _FOCUS)], darkcolor=[("focus", _FOCUS)])
    for field in ("TSpinbox", "TCombobox"):
        style.configure(field, background=_BUTTON, arrowsize=11, arrowcolor=MUTED)
        style.map(field, background=[("disabled", BACKGROUND), ("active", _BUTTON_HOVER)],
                  arrowcolor=[("disabled", DISABLED)])
    root.option_add("*TCombobox*Listbox.font", UI_FONT)
    root.option_add("*TCombobox*Listbox.background", PANEL)
    root.option_add("*TCombobox*Listbox.foreground", INK)
    root.option_add("*TCombobox*Listbox.selectBackground", SELECTION)
    root.option_add("*TCombobox*Listbox.selectForeground", INK)

    style.configure("TCheckbutton", background=BACKGROUND, foreground=INK, padding=(0, 2),
                    indicatorbackground=PANEL, indicatorforeground=ACCENT)
    style.configure("Panel.TCheckbutton", background=PANEL)
    style.map("TCheckbutton", indicatorbackground=[("disabled", HEADER), ("selected", SELECTION)],
              indicatorforeground=[("disabled", DISABLED), ("selected", ACCENT)])
    style.configure("TNotebook", background=BACKGROUND, borderwidth=1, tabmargins=(0, 0, 0, 0))
    style.configure("TNotebook.Tab", padding=(13, 5), background=HEADER, foreground=MUTED, borderwidth=1)
    style.map("TNotebook.Tab", background=[("selected", PANEL), ("active", _BUTTON_HOVER)],
              foreground=[("selected", ACCENT)], padding=[("selected", (13, 5)), ("!selected", (13, 5))])

    line_height = tkfont.Font(root=root, font=UI_FONT).metrics("linespace")
    style.configure("Treeview", background=PANEL, fieldbackground=PANEL, foreground=INK,
                    rowheight=line_height + 6, borderwidth=1, relief="flat")
    style.configure("Treeview.Heading", background=HEADER, foreground=INK, relief="flat",
                    padding=(6, 4), font=(UI_FONT[0], UI_FONT[1], "bold"))
    style.map("Treeview", background=[("selected", SELECTION)], foreground=[("selected", INK)])
    style.map("Treeview.Heading", background=[("pressed", _BUTTON_PRESSED), ("active", _BUTTON_HOVER)])
    style.configure("TPanedwindow", background=BACKGROUND, sashwidth=6, sashrelief="flat")
    style.configure("Sash", sashthickness=6)
    style.configure("TScrollbar", background=_THUMB, troughcolor=HEADER, arrowcolor=MUTED,
                    borderwidth=0, arrowsize=11, gripcount=0, relief="flat",
                    lightcolor=_THUMB, darkcolor=_THUMB)
    style.map("TScrollbar", background=[("active", AXIS)])
    style.configure("Horizontal.TProgressbar", background=ACCENT, troughcolor=HEADER,
                    lightcolor=ACCENT, darkcolor=ACCENT, thickness=6, borderwidth=1)


def style_review_axis(axis):
    axis.set_facecolor(PANEL)
    axis.spines[["top", "right"]].set_visible(False)
    axis.spines[["bottom", "left"]].set_color(AXIS)
    axis.spines[["bottom", "left"]].set_linewidth(.75)
    axis.set_axisbelow(True)
    axis.grid(axis="y", color=RULE, linewidth=.6)
    axis.tick_params(colors=MUTED, labelsize=9, length=3, width=.7)
