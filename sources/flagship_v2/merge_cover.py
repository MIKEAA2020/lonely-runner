#!/usr/bin/env python3
"""Merge cover.pdf (page 0) with main.pdf (body) -> final paper PDF,
scaling the cover page EXACTLY to the body's page size (A4 points)."""
from pypdf import PdfReader, PdfWriter

body = PdfReader('main.pdf')
target_w = float(body.pages[0].mediabox.width)
target_h = float(body.pages[0].mediabox.height)

cover = PdfReader('cover.pdf')
cover_page = cover.pages[0]
cw, ch = float(cover_page.mediabox.width), float(cover_page.mediabox.height)
if abs(cw - target_w) > 0.05 or abs(ch - target_h) > 0.05:
    cover_page.scale_to(target_w, target_h)

writer = PdfWriter()
writer.add_page(cover_page)
for page in body.pages:
    writer.add_page(page)
with open('final.pdf', 'wb') as f:
    writer.write(f)
print('merged -> final.pdf, pages:', len(writer.pages))

