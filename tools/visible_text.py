#!/usr/bin/env python3
"""Write the visible text of a manuscript .docx, one body paragraph per line, as ``<index>\\t<text>``.

This is how the reference text ``visible_v33.txt`` used by ``tools/check_numbers.py`` is made
(``visible_v30f.txt`` was made the same way). Every ``w:p`` element of the document body is one
line, numbered from 0 in document order (table-cell paragraphs included). Text inside tracked
deletions (``w:del``) is skipped, tracked insertions are kept, and tabs become spaces.

Usage
-----
    python tools/visible_text.py FastaLake_v33.docx visible_v33.txt
"""

import sys
import zipfile
import xml.etree.ElementTree as ET

W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"


def _text(el, in_del=False):
    """Concatenate w:t and w:tab descendants of ``el``, skipping those inside w:del."""
    out = []
    for ch in el:
        d = in_del or ch.tag == W + "del"
        if ch.tag == W + "t" and not d:
            out.append(ch.text or "")
        elif ch.tag == W + "tab" and not d:
            out.append(" ")
        out.append(_text(ch, d))
    return "".join(out)


def visible_paragraphs(docx):
    """Return the visible text of every body paragraph of ``docx`` in document order."""
    root = ET.fromstring(zipfile.ZipFile(docx).read("word/document.xml"))
    return [_text(p) for p in root.find(W + "body").iter(W + "p")]


def main(argv=None):
    """Command-line entry point."""
    a = argv if argv is not None else sys.argv[1:]
    if len(a) != 2:
        sys.exit(__doc__)
    with open(a[1], "w", encoding="utf-8") as f:
        for i, s in enumerate(visible_paragraphs(a[0])):
            f.write(f"{i}\t{s}\n")


if __name__ == "__main__":
    main()
