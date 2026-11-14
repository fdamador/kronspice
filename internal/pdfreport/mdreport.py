# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Render an example's Markdown (README.md or REPORT.md) into a kronreport.Report.

For examples whose report is the README itself, so the PDF and the README never
disagree. Used by examples/<name>/viz/make_report.py:

    from kronreport import Report
    from mdreport import render_markdown, add_figures, add_run_excerpt

    r = Report(...)
    used = render_markdown(r, os.path.join(EX, "README.md"), base=EX)
    add_figures(r, os.path.join(EX, "viz"), skip=used)
    add_run_excerpt(r, os.path.join(EX, "run_output.txt"))
    r.build(out)

Handled: headings (# is skipped as the title, ## → h1, ### and deeper → h2),
paragraphs, bullet and numbered lists (nested items flattened), pipe tables,
fenced code, images (PNG/JPG/GIF that exist), **bold**, *italic*, `code`,
[links](url) (text kept), and Unicode super/subscript digits turned into
<super>/<sub> markup so Arial renders them.
"""
import os
import re
from xml.sax.saxutils import escape as esc

__all__ = ["render_markdown", "add_figures", "add_run_excerpt", "inline"]

_SUP = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺", "0123456789−+")
_SUB = str.maketrans("₀₁₂₃₄₅₆₇₈₉", "0123456789")
# glyphs Arial lacks (kronreport's notes): spelled out
_SWAP = {"∝": " ∼ ", "∇": "grad ", "⟨": "<", "⟩": ">", "≲": "≤", "∥": "||", "⊥": " perp ", "∘": "°",
         "∓": "−/+", "∠": "angle "}


def _scripts(s):
    s = re.sub(r"[⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺]+", lambda m: "<super>" + m.group(0).translate(_SUP) + "</super>", s)
    s = re.sub(r"[₀₁₂₃₄₅₆₇₈₉]+", lambda m: "<sub>" + m.group(0).translate(_SUB) + "</sub>", s)
    return s


def inline(text):
    """Markdown inline → reportlab paragraph markup."""
    for k, v in _SWAP.items():
        text = text.replace(k, v)
    codes = []

    def keep(m):
        codes.append(m.group(1))
        return f"\x00{len(codes) - 1}\x00"

    text = re.sub(r"`([^`]+)`", keep, text)
    text = re.sub(r"!\[([^\]]*)\]\([^)]*\)", r"\1", text)
    text = re.sub(r"\[([^\]]+)\]\([^)]*\)", r"\1", text)
    text = esc(text)
    text = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", text)
    text = re.sub(r"(?<![\w*])\*(?!\s)(.+?)(?<!\s)\*(?![\w*])", r"<i>\1</i>", text)
    text = re.sub(r"(?<![\w_])_(?!\s)([^_]+?)(?<!\s)_(?![\w_])", r"<i>\1</i>", text)
    text = _scripts(text)
    text = re.sub(r"\x00(\d+)\x00", lambda m: "<i>" + esc(codes[int(m.group(1))]) + "</i>", text)
    return text


def _cells(line):
    line = line.strip()
    if line.startswith("|"):
        line = line[1:]
    if line.endswith("|"):
        line = line[:-1]
    return [c.strip() for c in re.split(r"(?<!\\)\|", line)]


def _widths(rows):
    n = max(len(r) for r in rows)
    w = [0.0] * n
    for r in rows:
        for i, c in enumerate(r):
            w[i] = max(w[i], min(len(re.sub(r"<[^>]+>", "", c)), 60))
    w = [max(x, 4) for x in w]
    t = sum(w)
    return [x / t for x in w]


def render_markdown(r, path, base=None, skip_title=True):
    """Render the file into r; returns the set of image paths it drew."""
    base = base or os.path.dirname(path)
    lines = open(path, encoding="utf-8").read().splitlines()
    used = set()
    para, items = [], []
    i = 0
    seen_title = False

    def flush():
        if para:
            r.p(inline(" ".join(para)))
            para.clear()
        if items:
            r.bullets([inline(t) for t in items])
            items.clear()

    while i < len(lines):
        ln = lines[i]
        s = ln.strip()
        if s.startswith("```"):
            flush()
            j = i + 1
            buf = []
            while j < len(lines) and not lines[j].strip().startswith("```"):
                buf.append(lines[j])
                j += 1
            r.code("\n".join(buf))
            i = j + 1
            continue
        m = re.match(r"^(#{1,6})\s+(.*)$", s)
        if m:
            flush()
            lvl, t = len(m.group(1)), inline(m.group(2))
            if lvl == 1 and skip_title and not seen_title:
                seen_title = True
            elif lvl <= 2:
                r.h1(t)
            else:
                r.h2(t)
            i += 1
            continue
        if s.startswith("|") and i + 1 < len(lines) and re.match(r"^\|?\s*:?-{2,}", lines[i + 1].strip()):
            flush()
            rows = [_cells(s)]
            j = i + 2
            while j < len(lines) and lines[j].strip().startswith("|"):
                rows.append(_cells(lines[j]))
                j += 1
            n = len(rows[0])
            rows = [(row + [""] * n)[:n] for row in rows]
            rows = [[inline(c) for c in row] for row in rows]
            r.table(rows, _widths(rows))
            i = j
            continue
        m = re.match(r"^!\[([^\]]*)\]\(([^)\s]+)", s)
        if m:
            flush()
            p = os.path.normpath(os.path.join(base, m.group(2)))
            if os.path.splitext(p)[1].lower() in (".png", ".jpg", ".jpeg", ".gif") and os.path.exists(p):
                r.figure(p, inline(m.group(1) or os.path.basename(p)), max_height_mm=150)
                used.add(p)
            i += 1
            continue
        m = re.match(r"^(\s*)(?:[-*+]|\d+[.)])\s+(.*)$", ln)
        if m and s:
            if para:
                r.p(inline(" ".join(para)))
                para.clear()
            items.append(m.group(2))
            i += 1
            continue
        if not s:
            flush()
            i += 1
            continue
        if items and (ln.startswith("  ") or ln.startswith("\t")):
            items[-1] += " " + s
        else:
            if items:
                r.bullets([inline(t) for t in items])
                items.clear()
            para.append(s)
        i += 1
    flush()
    return used


def add_figures(r, vizdir, skip=(), heading="Figures"):
    """Add every PNG/JPG (and GIF first frame) in vizdir not already drawn."""
    if not os.path.isdir(vizdir):
        return
    skip = {os.path.normpath(p) for p in skip}
    names = sorted(f for f in os.listdir(vizdir) if os.path.splitext(f)[1].lower() in (".png", ".jpg", ".jpeg", ".gif"))
    # a GIF whose PNG still exists is the same picture
    stems = {os.path.splitext(f)[0] for f in names if not f.lower().endswith(".gif")}
    names = [f for f in names if not (f.lower().endswith(".gif") and os.path.splitext(f)[0] in stems)]
    names = [f for f in names if os.path.normpath(os.path.join(vizdir, f)) not in skip]
    if not names:
        return
    r.h1(heading)
    for k, f in enumerate(names, 1):
        cap = f"Figure {k}. {esc(os.path.splitext(f)[0].replace('_', ' '))} (viz/{esc(f)})"
        r.figure(os.path.join(vizdir, f), cap, max_height_mm=150)


def add_run_excerpt(r, path, max_lines=120, heading="Run output (excerpt)"):
    """The head of run_output.txt, as recorded (not re-run)."""
    if not os.path.exists(path):
        return
    lines = open(path, encoding="utf-8", errors="replace").read().splitlines()
    r.h1(heading)
    note = f"First {min(max_lines, len(lines))} of {len(lines)} lines of {os.path.basename(path)}, as recorded."
    r.small(esc(note))
    text = "\n".join(lines[:max_lines])
    for k, v in _SWAP.items():
        text = text.replace(k, v)
    r.code(text)
