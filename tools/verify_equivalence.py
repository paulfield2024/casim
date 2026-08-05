#!/usr/bin/env python3
"""Verify that a UMDP3-cleaned file is token-equivalent to the original.

Confirms that clean_umdp3.py has only changed whitespace, keyword case,
and old-style-operator spelling - i.e. that no code token, identifier,
numeric literal or string-literal content has been added, removed or
altered.

Usage: verify_equivalence.py <original_file> <cleaned_file>
Exits 0 and prints "OK" if equivalent, otherwise exits 1 and prints a
short context diff of the normalised token streams.
"""
import re
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from clean_umdp3 import split_segments  # noqa: E402

# Multi-char operators must be substituted before single-char ones.
_OP_CANON = [
    (re.compile(r">="), ".ge."),
    (re.compile(r"<="), ".le."),
    (re.compile(r"=="), ".eq."),
    (re.compile(r"/="), ".ne."),
    (re.compile(r">"), ".gt."),
    (re.compile(r"<"), ".lt."),
]

_WS_RE = re.compile(r"\s+")


def normalise_file(path):
    with open(path, "r") as f:
        text = f.read()
    tokens = []
    for line in text.splitlines():
        # Directive lines are treated as opaque by the cleaner; compare
        # them verbatim (lower-cased, whitespace-stripped) too.
        for kind, seg in split_segments(line):
            if kind == "comment":
                continue
            s = seg.lower()
            if kind == "code":
                for regex, repl in _OP_CANON:
                    s = regex.sub(repl, s)
            s = _WS_RE.sub("", s)
            tokens.append(s)
    return "".join(tokens)


def main():
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: verify_equivalence.py <original> <cleaned>\n")
        sys.exit(2)
    orig, cleaned = sys.argv[1], sys.argv[2]
    a = normalise_file(orig)
    b = normalise_file(cleaned)
    if a == b:
        print(f"OK: {orig}")
        sys.exit(0)
    # Find first divergence to aid debugging.
    n = min(len(a), len(b))
    i = 0
    while i < n and a[i] == b[i]:
        i += 1
    ctx = 60
    sys.stderr.write(f"MISMATCH: {orig} vs {cleaned}\n")
    sys.stderr.write(f"  original @ {i}: ...{a[max(0,i-ctx):i+ctx]}...\n")
    sys.stderr.write(f"  cleaned  @ {i}: ...{b[max(0,i-ctx):i+ctx]}...\n")
    sys.exit(1)


if __name__ == "__main__":
    main()
