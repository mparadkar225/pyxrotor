"""List how each Fortran MIN/MAX in the compiled reference treats a NaN argument.

gfortran leaves MIN/MAX of a NaN unspecified. At -O1 on x86-64 each call compiles to
MINSD/MAXSD (or a compare-and-select sequence), and those return their *source* operand
when either operand is NaN. Which Fortran argument ends up as the source differs from
site to site, so this tool recompiles the given sources with -g (checking the code is
byte-identical to the build without it), disassembles them, and prints every site with
its source line and operands. A memory source operand is a constant (decoded and shown):
that operand (a constant, a COMMON variable or a stack value) is the result whenever either
argument is NaN. Register-register sites
need reading the surrounding code (use --context); pin them with a non-finite test.

    python tools/minmax_sites.py xoper xdesi          # routines in those files
    python tools/minmax_sites.py xoper --routine APINIT --context 6
"""
from __future__ import annotations

import argparse
import re
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "reference" / "xrotor-7.55" / "src"
SHIM = ROOT / "reference" / "shim"
FLAGS = ["-O1", "-fPIC", "-fdefault-real-8", "-fdefault-double-8", "-ffp-contract=off",
         "-std=legacy", "-fallow-argument-mismatch", "-w"]


def compile_pair(path: Path, out: Path) -> Path:
    g = out / (path.stem + ".g.o")
    plain = out / (path.stem + ".o")
    for obj, extra in ((g, ["-g"]), (plain, [])):
        subprocess.run(["gfortran", *FLAGS, *extra, f"-I{SRC}", "-c", str(path), "-o", str(obj)], check=True)
    dis = [subprocess.run(["objdump", "-d", "--no-show-raw-insn", str(o)], capture_output=True,
                          text=True, check=True).stdout.splitlines()[3:] for o in (g, plain)]
    if dis[0] != dis[1]:
        raise SystemExit(f"{path.name}: -g changed the generated code; results would not apply")
    return g


def sites(obj: Path, routine: str | None, context: int):
    dis = subprocess.run(["objdump", "-dlr", "--no-show-raw-insn", str(obj)], capture_output=True,
                         text=True, check=True).stdout.splitlines()
    func, line = None, None
    lines_of = {}
    for k, text in enumerate(dis):
        m = re.match(r"^[0-9a-f]+ <(\w+)>:", text)
        if m:
            func = m.group(1)
            continue
        m = re.match(r"^(/\S+):(\d+)", text)
        if m:
            line = (m.group(1), int(m.group(2)))
            continue
        if func is None or (routine and func.rstrip("_").upper() != routine.upper()):
            continue
        m = re.search(r"\b(maxsd|minsd)\s+(\S+),(%xmm\d+)", text)
        if not m:
            continue
        op, srcop, dst = m.groups()
        yield func, line, op, srcop, dst, dis[max(0, k - context):k + 1]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+", help="source files without .f (xoper, xdesi, ...; 'shim' for the shim)")
    ap.add_argument("--routine")
    ap.add_argument("--context", type=int, default=0)
    a = ap.parse_args()
    with tempfile.TemporaryDirectory() as tmp:
        for name in a.files:
            path = SHIM / "ref755_shim.f" if name == "shim" else SRC / f"{name}.f"
            obj = compile_pair(path, Path(tmp))
            src = path.read_text(errors="replace").splitlines()
            for func, line, op, srcop, dst, ctx in sites(obj, a.routine, a.context):
                ftext = src[line[1] - 1].strip() if line and Path(line[0]).name == path.name else "?"
                ln = line[1] if line else "?"
                nan = ("NaN -> constant" if "(%rip)" in srcop else
                       "NaN -> memory operand" if "(" in srcop else f"NaN -> source register {srcop}")
                print(f"{func.rstrip('_').upper():10s} {path.name}:{ln:<5} {op} {srcop},{dst}  [{nan}]  | {ftext}")
                if a.context:
                    for c in ctx:
                        print("      " + c)


if __name__ == "__main__":
    main()
