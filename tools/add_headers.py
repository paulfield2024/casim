#!/usr/bin/env python3
"""Prepend a standard UMDP3-style header block to each cleaned file.

Purely additive: inserts a fixed comment block (copyright + one-line
description + code-owner/language boilerplate) at the very top of the
file. No existing line is altered, moved or removed.

Usage: add_headers.py <clean_src_root> <descriptions.json>
Rewrites every .F90 file under <clean_src_root> in place.
"""
import json
import sys
from pathlib import Path

HEADER_TEMPLATE = """\
! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   {description}
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
"""


def main():
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: add_headers.py <clean_src_root> <descriptions.json>\n")
        sys.exit(2)
    root = Path(sys.argv[1])
    with open(sys.argv[2]) as f:
        descriptions = json.load(f)

    missing = []
    for f in sorted(root.rglob("*.F90")):
        rel = str(f.relative_to(root))
        desc = descriptions.get(rel)
        if desc is None:
            missing.append(rel)
            continue
        text = f.read_text()
        header = HEADER_TEMPLATE.format(description=desc)
        f.write_text(header + text)

    if missing:
        sys.stderr.write("No description for:\n  " + "\n  ".join(missing) + "\n")
        sys.exit(1)
    print(f"Headers added to files under {root}")


if __name__ == "__main__":
    main()
