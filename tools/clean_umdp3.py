#!/usr/bin/env python3
"""UMDP3 mechanical style cleaner for CASIM Fortran source.

Applies a small set of *safe*, purely-lexical UMDP3 (section 3.4/3.9)
formatting rules to each .F90 file:

  - old-style relational operators (.EQ. .NE. .GT. .LT. .GE. .LE.)
    replaced with their symbolic equivalents (== /= > < >= <=)
  - logical operators/literals (.AND. .OR. .NOT. .EQV. .NEQV.
    .TRUE. .FALSE.) upper-cased in dot-form
  - compressed keyword forms (ENDIF, ENDDO, ELSEIF, ELSEWHERE,
    ENDSELECT, ENDWHERE, ENDTYPE, ENDMODULE, ENDSUBROUTINE,
    ENDFUNCTION, ENDINTERFACE, ENDPROGRAM, ENDASSOCIATE, ENDBLOCK,
    GOTO, DOUBLEPRECISION) expanded to the two-word spaced form
  - a curated list of Fortran structural keywords upper-cased
    wherever they appear as whole words in executable/declaration
    code (never inside string literals or comments, and never
    immediately preceded by '%' so derived-type components are
    left untouched)

Only whitespace/case is changed by this script - no code tokens,
identifiers, literals (other than the operator forms above), or
comments are altered, added or removed except where a new standard
header is prepended (handled separately by add_headers.py).

Usage: clean_umdp3.py <src_file> <dst_file>
"""
import re
import sys

# Compressed/joined keyword forms -> spaced, upper-case two-word forms.
# Matched as whole words (case-insensitive) so e.g. "myenddo" is untouched.
COMPRESSED = [
    (r"endif", "END IF"),
    (r"enddo", "END DO"),
    (r"elseif", "ELSE IF"),
    (r"elsewhere", "ELSE WHERE"),
    (r"endselect", "END SELECT"),
    (r"endwhere", "END WHERE"),
    (r"endtype", "END TYPE"),
    (r"endmodule", "END MODULE"),
    (r"endsubroutine", "END SUBROUTINE"),
    (r"endfunction", "END FUNCTION"),
    (r"endinterface", "END INTERFACE"),
    (r"endprogram", "END PROGRAM"),
    (r"endassociate", "END ASSOCIATE"),
    (r"endblock", "END BLOCK"),
    (r"goto", "GO TO"),
    (r"doubleprecision", "DOUBLE PRECISION"),
]

# Old-style relational operators -> symbolic form (UMDP3 3.9).
RELATIONAL = [
    (r"\.eq\.", "=="),
    (r"\.ne\.", "/="),
    (r"\.ge\.", ">="),
    (r"\.le\.", "<="),
    (r"\.gt\.", ">"),
    (r"\.lt\.", "<"),
]

# Logical operators/literals -> upper-cased dot-form (UMDP3 3.4).
LOGICAL_DOTS = [
    (r"\.and\.", ".AND."),
    (r"\.or\.", ".OR."),
    (r"\.not\.", ".NOT."),
    (r"\.eqv\.", ".EQV."),
    (r"\.neqv\.", ".NEQV."),
    (r"\.true\.", ".TRUE."),
    (r"\.false\.", ".FALSE."),
]

# Structural Fortran keywords to upper-case wherever they appear as a
# whole word in code (not comments/strings). Deliberately excludes
# generic intrinsic-function-only names and anything not confirmed
# safe against this codebase's identifier usage.
KEYWORDS = [
    "module", "subroutine", "function", "program", "interface",
    "contains", "implicit", "none", "private", "public", "parameter",
    "allocatable", "pointer", "target", "optional", "save", "dimension",
    "intent", "kind", "type", "class", "real", "integer", "logical",
    "character", "complex", "double", "precision", "use", "only",
    "call", "return", "if", "then", "else", "do", "while", "cycle",
    "exit", "select", "case", "where", "allocate", "deallocate",
    "associate", "block", "enum", "procedure", "operator", "assignment",
    "continue", "stop", "format", "write", "read", "open", "close",
    "print", "namelist", "common", "equivalence", "external",
    "intrinsic", "data", "end", "in", "out", "inout", "recursive",
    "elemental", "pure", "result", "entry", "go", "to",
]
KEYWORD_RE = re.compile(
    r"(?<!%)\b(" + "|".join(KEYWORDS) + r")\b", re.IGNORECASE
)
COMPRESSED_RES = [
    (re.compile(r"\b" + pat + r"\b", re.IGNORECASE), repl)
    for pat, repl in COMPRESSED
]
RELATIONAL_RES = [(re.compile(pat, re.IGNORECASE), repl) for pat, repl in RELATIONAL]
LOGICAL_RES = [(re.compile(pat, re.IGNORECASE), repl) for pat, repl in LOGICAL_DOTS]


def split_segments(line):
    """Split a line into (kind, text) segments, kind in
    {'code', 'string', 'comment'}.

    'string' segments (including their delimiting quotes) and the
    trailing 'comment' segment (starting at the first unquoted '!')
    are returned verbatim so they are never touched by later
    transforms; only 'code' segments are transformable.
    """
    segments = []
    buf = []
    i = 0
    n = len(line)
    while i < n:
        ch = line[i]
        if ch in ("'", '"'):
            if buf:
                segments.append(("code", "".join(buf)))
                buf = []
            quote = ch
            j = i + 1
            while j < n:
                if line[j] == quote:
                    # Doubled quote is an escaped literal quote.
                    if j + 1 < n and line[j + 1] == quote:
                        j += 2
                        continue
                    j += 1
                    break
                j += 1
            segments.append(("string", line[i:j]))
            i = j
            continue
        if ch == "!":
            if buf:
                segments.append(("code", "".join(buf)))
                buf = []
            segments.append(("comment", line[i:]))
            return segments
        buf.append(ch)
        i += 1
    if buf:
        segments.append(("code", "".join(buf)))
    return segments


def clean_code(code):
    for regex, repl in COMPRESSED_RES:
        code = regex.sub(repl, code)
    for regex, repl in RELATIONAL_RES:
        code = regex.sub(repl, code)
    for regex, repl in LOGICAL_RES:
        code = regex.sub(repl, code)
    code = KEYWORD_RE.sub(lambda m: m.group(1).upper(), code)
    return code


def clean_line(line):
    # Preserve the line ending exactly.
    ending = ""
    if line.endswith("\r\n"):
        line, ending = line[:-2], "\r\n"
    elif line.endswith("\n"):
        line, ending = line[:-1], "\n"

    line = line.expandtabs(8)

    # A directive line (e.g. "!$OMP ...") is a full-line comment as far
    # as a non-OpenMP-aware transform is concerned; leave untouched.
    # C-preprocessor directives ("#if", "#else", "#endif", ...) are
    # case-sensitive to cpp, not Fortran keywords; leave untouched too.
    if re.match(r"^\s*!\$", line) or re.match(r"^\s*#", line):
        return line + ending

    out = []
    for kind, text in split_segments(line):
        out.append(clean_code(text) if kind == "code" else text)
    return "".join(out) + ending


def clean_text(text):
    lines = text.splitlines(keepends=True)
    return "".join(clean_line(l) for l in lines)


def main():
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: clean_umdp3.py <src_file> <dst_file>\n")
        sys.exit(2)
    src, dst = sys.argv[1], sys.argv[2]
    with open(src, "r") as f:
        text = f.read()
    cleaned = clean_text(text)
    with open(dst, "w") as f:
        f.write(cleaned)


if __name__ == "__main__":
    main()
