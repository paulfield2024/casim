#!/usr/bin/env python3
"""Add CASIM-paper documentation (module-level Method/Paper reference
block + per-routine citation comments) to every file under a cleaned
source tree, writing the result to a new tree. Purely additive: only
inserts new comment lines, never edits, removes or reorders any existing
line, so the resulting files are token-equivalent to the input (see
verify_equivalence.py) and produce identical compiled behaviour.

Usage: add_paper_docs.py <src_clean_root> <src_clean_documented_root>
"""
import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paper_docs import DOCS  # noqa: E402

HEADER_ANCHOR = (
    "!   This code follows the Met Office UMDP3 style guide.\n!\n"
)

ROUTINE_RE = re.compile(
    r"^(?P<indent>[ \t]*)(?P<prefix>(?:RECURSIVE\s+)?)"
    r"(?P<kind>SUBROUTINE|FUNCTION)\s+(?P<name>\w+)",
    re.IGNORECASE,
)


def wrap_block(title, lines):
    out = [f"! {title}:\n"]
    for line in lines:
        out.append(f"!   {line}\n")
    out.append("!\n")
    return out


def insert_method_block(text, doc):
    if HEADER_ANCHOR not in text:
        raise ValueError("header anchor not found - was add_headers.py run first?")
    block = "".join(wrap_block("Method", doc["method"]) + wrap_block(
        "Paper reference", doc["reference"]))
    return text.replace(HEADER_ANCHOR, HEADER_ANCHOR + block, 1)


def insert_routine_comments(text, routines):
    if not routines:
        return text
    out_lines = []
    seen = set()
    for line in text.splitlines(keepends=True):
        m = ROUTINE_RE.match(line)
        if m:
            name = m.group("name").upper()
            doc_lines = routines.get(name)
            if doc_lines and name not in seen:
                seen.add(name)
                indent = m.group("indent")
                for dl in doc_lines:
                    out_lines.append(f"{indent}! {dl}\n")
        out_lines.append(line)
    return "".join(out_lines)


def main():
    if len(sys.argv) != 3:
        sys.stderr.write(
            "Usage: add_paper_docs.py <src_clean_root> <src_clean_documented_root>\n"
        )
        sys.exit(2)
    src_root = Path(sys.argv[1])
    dst_root = Path(sys.argv[2])

    if dst_root.exists():
        shutil.rmtree(dst_root)
    shutil.copytree(src_root, dst_root)

    missing = []
    for f in sorted(dst_root.rglob("*.F90")):
        rel = str(f.relative_to(dst_root))
        doc = DOCS.get(rel)
        if doc is None:
            missing.append(rel)
            continue
        text = f.read_text()
        text = insert_method_block(text, doc)
        text = insert_routine_comments(text, doc.get("routines", {}))
        f.write_text(text)

    if missing:
        sys.stderr.write("No paper-doc entry for:\n  " + "\n  ".join(missing) + "\n")
        sys.exit(1)
    print(f"Paper documentation added to files under {dst_root}")


if __name__ == "__main__":
    main()
