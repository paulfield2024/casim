#!/usr/bin/env python3
"""Topologically sort CASIM Fortran source files by module dependency.

Scans a set of source directories for .F90/.f90 files, extracts the module(s)
each file defines (`module X`) and the modules each file uses (`use Y`), then
prints a build order (one path per line) such that every file appears after
all the files defining modules it uses.

Usage: build_order.py dir1 [dir2 ...]
"""
import re
import sys
from pathlib import Path

MODULE_RE = re.compile(r"^\s*module\s+(\w+)\s*$", re.IGNORECASE)
USE_RE = re.compile(r"^\s*use\s*(?:,\s*intrinsic\s*)?(?:::)?\s*(\w+)", re.IGNORECASE)


def scan(paths):
    files = []
    for p in paths:
        d = Path(p)
        files.extend(sorted(d.rglob("*.F90")))
        files.extend(sorted(d.rglob("*.f90")))
    return files


def parse(files):
    provides = {}   # module_name -> file
    needs = {}      # file -> set(module_names)
    for f in files:
        mods = set()
        text = f.read_text(errors="ignore")
        for line in text.splitlines():
            m = MODULE_RE.match(line)
            if m and m.group(1).lower() != "procedure":
                provides[m.group(1).lower()] = f
            u = USE_RE.match(line)
            if u:
                mods.add(u.group(1).lower())
        needs[f] = mods
    return provides, needs


def topo_sort(files, provides, needs):
    order = []
    visited = set()
    visiting = set()

    def visit(f):
        if f in visited:
            return
        if f in visiting:
            raise RuntimeError(f"Circular dependency involving {f}")
        visiting.add(f)
        for mod in needs.get(f, ()):
            dep_file = provides.get(mod)
            if dep_file is not None and dep_file != f:
                visit(dep_file)
        visiting.discard(f)
        visited.add(f)
        order.append(f)

    for f in files:
        visit(f)
    return order


if __name__ == "__main__":
    dirs = sys.argv[1:]
    if not dirs:
        print(__doc__)
        sys.exit(1)
    files = scan(dirs)
    provides, needs = parse(files)
    order = topo_sort(files, provides, needs)
    for f in order:
        print(f)
