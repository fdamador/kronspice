# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Shared layout for the examples' PDF reports (reportlab, A4).

Each example has examples/<name>/viz/make_report.py, which imports this module:

    import os, sys
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    "..", "..", "..", "internal", "pdfreport"))
    from kronreport import *

    r = Report("examples/<name>", "Title", "subtitle", "register entry / CLAUDE.md section")
    r.h1("Summary"); r.p("..."); r.bullets([...])
    r.table([["head", ...], [...]], [0.4, 0.2, ...], bold_rows=(2,))
    r.figure(os.path.join(HERE, "x.png"), "Figure 1. ...")
    r.code("go run ./examples/<name>")
    r.build(out_path)

Text is reportlab paragraph markup: <b>, <i>, <sub>, <super>; escape &, <, >
in plain text with esc(). Use a Unicode TTF (Arial on Windows, DejaVu Sans
elsewhere) so Greek letters, ×, ≤, Ω, µ, ² and − render. Write exponents as
10<super>−8</super>, never with Unicode superscript characters such as ⁻⁸.
Arial has no glyph for ∇ ∝ ∠ ∘ ∓ ∥ ⊥ ≲ ⟨ ⟩ or subscript digits (₀): write them
in words or markup. The code() font (Consolas) has no Greek.
"""
import os
import sys
from datetime import date
from xml.sax.saxutils import escape as esc

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (Image, KeepTogether, PageBreak, Paragraph, Preformatted,
                                SimpleDocTemplate, Spacer, Table, TableStyle)

__all__ = ["Report", "esc", "mm", "COPYRIGHT"]

COPYRIGHT = "© 2026 fdamador@comcast.net"

_FONTS = [
    ("C:/Windows/Fonts/arial.ttf", "C:/Windows/Fonts/arialbd.ttf", "C:/Windows/Fonts/ariali.ttf"),
    ("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
     "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
     "/usr/share/fonts/truetype/dejavu/DejaVuSans-Oblique.ttf"),
    ("/Library/Fonts/Arial.ttf", "/Library/Fonts/Arial Bold.ttf", "/Library/Fonts/Arial Italic.ttf"),
]
_MONO = ["C:/Windows/Fonts/consola.ttf", "/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf"]


def _register():
    if "U" in pdfmetrics.getRegisteredFontNames():
        return
    for r, b, i in _FONTS:
        if os.path.exists(r):
            pdfmetrics.registerFont(TTFont("U", r))
            pdfmetrics.registerFont(TTFont("U-B", b))
            pdfmetrics.registerFont(TTFont("U-I", i))
            pdfmetrics.registerFontFamily("U", normal="U", bold="U-B", italic="U-I", boldItalic="U-B")
            break
    else:
        sys.exit("kronreport: no Unicode TTF font found (Arial or DejaVu Sans)")
    for m in _MONO:
        if os.path.exists(m):
            pdfmetrics.registerFont(TTFont("UM", m))
            return
    raise SystemExit("kronreport: no Unicode monospace font found (Consolas or DejaVu Sans Mono)")


_register()
_ss = getSampleStyleSheet()
NAVY = colors.HexColor("#1f3a5f")
BODY = ParagraphStyle("body", parent=_ss["Normal"], fontName="U", fontSize=9.5, leading=13.2, spaceAfter=5)
SMALL = ParagraphStyle("small", parent=BODY, fontSize=8, leading=10.5, textColor=colors.HexColor("#444444"))
CAP = ParagraphStyle("cap", parent=SMALL, alignment=TA_CENTER, spaceBefore=3, spaceAfter=10)
H1 = ParagraphStyle("h1", parent=_ss["Heading1"], fontName="U-B", fontSize=14, leading=18,
                    spaceBefore=10, spaceAfter=6, textColor=NAVY, keepWithNext=1)
H2 = ParagraphStyle("h2", parent=_ss["Heading2"], fontName="U-B", fontSize=11, leading=14,
                    spaceBefore=8, spaceAfter=4, textColor=NAVY, keepWithNext=1)
TITLE = ParagraphStyle("title", parent=_ss["Title"], fontName="U-B", fontSize=19, leading=24, spaceAfter=4)
SUB = ParagraphStyle("sub", parent=BODY, alignment=TA_CENTER, fontSize=10.5, textColor=colors.HexColor("#555555"))
CELL = ParagraphStyle("cell", parent=BODY, fontSize=8, leading=10, spaceAfter=0)
CELLB = ParagraphStyle("cellb", parent=CELL, fontName="U-B")
BULLET = ParagraphStyle("bullet", parent=BODY, leftIndent=12, bulletIndent=2, spaceAfter=3)
CODE = ParagraphStyle("code", parent=BODY, fontName="UM", fontSize=7.6, leading=9.6,
                      backColor=colors.HexColor("#f3f4f6"), borderPadding=4, spaceBefore=4, spaceAfter=8)


class Report:
    """A report under construction. Widths are fractions of the text width."""

    def __init__(self, example, title, subtitle, register, run_note=""):
        self.example = example          # e.g. "examples/al_cell_pump"
        self.title = title
        self.register = register        # e.g. "register entry H38" or "CLAUDE.md §3e"
        self.W = A4[0] - 36 * mm
        self.s = [Spacer(1, 6 * mm), Paragraph(title, TITLE), Paragraph(subtitle, SUB)]
        line = f"kronspice simulation report · {date.today().isoformat()}"
        if run_note:
            line += " · " + run_note
        self.s += [Paragraph(line, SUB), Paragraph(COPYRIGHT, SUB), Spacer(1, 6 * mm)]

    # text
    def h1(self, t):
        self.s.append(Paragraph(t, H1))

    def h2(self, t):
        self.s.append(Paragraph(t, H2))

    def p(self, t):
        self.s.append(Paragraph(t, BODY))

    def small(self, t):
        self.s.append(Paragraph(t, SMALL))

    def bullets(self, items):
        self.s += [Paragraph(t, BULLET, bulletText="•") for t in items]

    def code(self, text):
        """Monospace block, e.g. run commands or a short excerpt of run_output.txt (plain text, not markup)."""
        self.s.append(Preformatted(text.rstrip("\n"), CODE))

    def space(self, h=3):
        self.s.append(Spacer(1, h * mm))

    def page_break(self):
        self.s.append(PageBreak())

    # tables and figures
    def table(self, rows, widths, bold_rows=()):
        """rows[0] is the header; cells are markup strings; widths are fractions summing to ~1."""
        data = [[Paragraph(str(c), CELLB if (r == 0 or r in bold_rows) else CELL) for c in row]
                for r, row in enumerate(rows)]
        t = Table(data, colWidths=[w * self.W for w in widths], repeatRows=1, hAlign="LEFT")
        st = [("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#dfe7f1")),
              ("LINEBELOW", (0, 0), (-1, 0), 0.6, NAVY),
              ("LINEBELOW", (0, -1), (-1, -1), 0.4, colors.HexColor("#999999")),
              ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
              ("TOPPADDING", (0, 0), (-1, -1), 2.5), ("BOTTOMPADDING", (0, 0), (-1, -1), 2.5)]
        for r in range(2, len(rows), 2):
            st.append(("BACKGROUND", (0, r), (-1, r), colors.HexColor("#f5f7fa")))
        for r in bold_rows:
            st.append(("BACKGROUND", (0, r), (-1, r), colors.HexColor("#fff4d6")))
        t.setStyle(TableStyle(st))
        self.s += [t, Spacer(1, 3 * mm)]

    def figure(self, path, caption, width_frac=1.0, max_height_mm=200):
        """PNG or JPG; a GIF is drawn from its first frame."""
        img = Image(path)
        w = width_frac * self.W
        h = w * img.imageHeight / img.imageWidth
        if h > max_height_mm * mm:
            h = max_height_mm * mm
            w = h * img.imageWidth / img.imageHeight
        img.drawWidth, img.drawHeight = w, h
        self.s.append(KeepTogether([img, Paragraph(caption, CAP)]))

    # output
    def _footer(self, canvas, doc):
        canvas.saveState()
        canvas.setFont("U", 7.5)
        canvas.setFillColor(colors.HexColor("#777777"))
        canvas.drawString(18 * mm, 10 * mm, f"kronspice · {self.example} · {self.register}")
        canvas.drawCentredString(A4[0] / 2, 6 * mm, COPYRIGHT)
        canvas.drawRightString(A4[0] - 18 * mm, 10 * mm, f"page {doc.page}")
        canvas.restoreState()

    def build(self, out):
        doc = SimpleDocTemplate(out, pagesize=A4, leftMargin=18 * mm, rightMargin=18 * mm,
                                topMargin=16 * mm, bottomMargin=18 * mm,
                                title=self.title, author="fdamador@comcast.net", subject=self.register)
        doc.build(self.s, onFirstPage=self._footer, onLaterPages=self._footer)
        print("wrote", out)
