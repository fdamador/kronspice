# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Stamp the copyright line onto an existing PDF that was not built with kronreport.

python internal/pdfreport/stamp_copyright.py FILE.pdf [FILE.pdf ...]

Each page gets the line kronreport.COPYRIGHT centred 6 mm above the bottom edge,
in 7.5 pt grey (the same footer the generated reports carry), and the document
metadata gets /Author and /Copyright. Pages already carrying the line are left
alone, so running it twice changes nothing. The file is rewritten in place.
"""
import io
import os
import sys

from pypdf import PdfReader, PdfWriter
from reportlab.lib import colors
from reportlab.lib.units import mm
from reportlab.pdfgen import canvas

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kronreport import COPYRIGHT  # noqa: E402  (also registers the "U" font)

EMAIL = "fdamador@comcast.net"


def overlay(w, h):
    buf = io.BytesIO()
    c = canvas.Canvas(buf, pagesize=(w, h))
    c.setFont("U", 7.5)
    c.setFillColor(colors.HexColor("#777777"))
    c.drawCentredString(w / 2, 6 * mm, COPYRIGHT)
    c.save()
    buf.seek(0)
    return PdfReader(buf).pages[0]


def stamp(path):
    r = PdfReader(path)
    w = PdfWriter(clone_from=r)
    n = 0
    for page in w.pages:
        if EMAIL in (page.extract_text() or ""):
            continue
        box = page.mediabox
        page.merge_page(overlay(float(box.width), float(box.height)))
        n += 1
    meta = dict(r.metadata or {})
    meta.update({"/Author": EMAIL, "/Copyright": COPYRIGHT})
    w.add_metadata(meta)
    with open(path, "wb") as f:
        w.write(f)
    print(f"{path}: stamped {n} of {len(w.pages)} pages")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    for p in sys.argv[1:]:
        stamp(p)
