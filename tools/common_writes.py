"""List the COMMON variables each XROTOR routine assigns to.

    python tools/common_writes.py SETXW TPQ APITER

Routines that INCLUDE 'XROTOR.INC' share every COMMON name, so a "local-looking" variable
such as RLX in SETXW is really solver state. Use this before porting a routine, and make
the port write every variable listed (or explain why not).
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from reference.commons import parse_include  # noqa: E402

SRC = ROOT / "reference" / "xrotor-7.55" / "src"
NAMES = {v.name for b in parse_include(SRC / "XROTOR.INC").values() for v in b.vars}


def routine_text(name: str) -> str:
    for f in sorted(SRC.glob("*.f")):
        lines = f.read_text(errors="replace").splitlines()
        for k, line in enumerate(lines):
            if re.match(rf"^      (SUBROUTINE|.*FUNCTION)\s+{name}\b", line, re.I):
                body = []
                for line2 in lines[k:]:
                    body.append(line2)
                    if re.match(r"^      END\s*(!.*)?$", line2) or re.match(r"^      END\s*$", line2):
                        break
                return "\n".join(body)
    raise SystemExit(f"{name} not found")


def common_writes(name: str) -> list[str]:
    text = routine_text(name)
    if "XROTOR.INC" not in text:
        return []
    writes = set()
    for line in text.splitlines():
        if not line.strip() or line[:1] in "Cc*!" or len(line) < 7:
            continue
        stmt = line[6:72]
        m = re.match(r"\s*(?:IF\s*\(.*\)\s*)?([A-Z][A-Z0-9_]*)\s*(\(.*?\))?\s*=(?!=)", stmt, re.I)
        if m and m.group(1).upper() in NAMES:
            writes.add(m.group(1).upper())
        m = re.match(r"\s*DO\s+\d*\s*([A-Z][A-Z0-9_]*)\s*=", stmt, re.I)   # DO-loop variables
        if m and m.group(1).upper() in NAMES:
            writes.add(m.group(1).upper())
    return sorted(writes)


if __name__ == "__main__":
    for r in sys.argv[1:]:
        print(f"{r}: {' '.join(common_writes(r.upper()))}")
