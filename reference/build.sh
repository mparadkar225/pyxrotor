#!/usr/bin/env bash
# Build the two XROTOR reference libraries into reference/build/.
#
#   reference/build/libxrotor755.so   XROTOR 7.55 (MIT source) + non-interactive shim
#   reference/build/xrotor-python/    DARcorporation xrotor-python, patched to double-precision I/O
#
# Requirements: gfortran, gcc, cmake, libX11 headers (the plot library is linked but never opened).
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC755="$HERE/xrotor-7.55"
SRCXP="$HERE/xrotor-python"
OUT="$HERE/build"

# Shared numerics for both references. Every REAL and DOUBLE PRECISION is 8 bytes, and fused
# multiply-add contraction is off so results are identical on x86-64 and ARM.
FFLAGS="-O1 -fPIC -fdefault-real-8 -fdefault-double-8 -ffp-contract=off"
F77FLAGS="$FFLAGS -std=legacy -fallow-argument-mismatch -w"

for d in "$SRC755/src" "$SRCXP/src"; do
  [ -d "$d" ] || { echo "Missing $d. Run: git submodule update --init" >&2; exit 1; }
done

rm -rf "$OUT"
mkdir -p "$OUT/obj755/plotlib" "$OUT/obj755/src" "$OUT/xrotor-python"

# ---------------------------------------------------------------- XROTOR 7.55
echo "== XROTOR 7.55: plot library"
( cd "$SRC755/plotlib"
  for f in plt_base plt_font plt_util plt_color set_subs gw_subs ps_subs plt_old plt_3D; do
    gfortran $F77FLAGS -c "$f.f" -o "$OUT/obj755/plotlib/$f.o"
  done
  gcc -O1 -fPIC -DUNDERSCORE -w -c Xwin.c -o "$OUT/obj755/plotlib/Xwin.o" )

echo "== XROTOR 7.55: solver"
# The main program becomes an ordinary subroutine so the library has no main().
sed 's/^      PROGRAM XROTOR/      SUBROUTINE XROTOR_MAIN/' "$SRC755/src/xrotor.f" > "$OUT/obj755/xrotor_lib.f"
SRCS="xoper xdesi xmodi xaero xjmap xio xnoise xrotpl xcasepl xbend xinte xutils jputil plutil
      modify srclin spline userio vortex plotdata"
gfortran $F77FLAGS -I"$SRC755/src" -c "$OUT/obj755/xrotor_lib.f" -o "$OUT/obj755/src/xrotor.o"
for f in $SRCS; do
  gfortran $F77FLAGS -I"$SRC755/src" -c "$SRC755/src/$f.f" -o "$OUT/obj755/src/$f.o"
done
gfortran $F77FLAGS -I"$SRC755/src" -c "$HERE/shim/ref755_shim.f" -o "$OUT/obj755/src/ref755_shim.o"

gfortran -shared -o "$OUT/libxrotor755.so" "$OUT"/obj755/src/*.o "$OUT"/obj755/plotlib/*.o -lX11

# ---------------------------------------------------------------- xrotor-python
echo "== xrotor-python"
cp -r "$SRCXP/src" "$SRCXP/xrotor" "$SRCXP/CMakeLists.txt" "$OUT/xrotor-python/"
# Its C interface passes single-precision floats. Switch the interface to double so
# comparisons are not limited to ~1e-7 by float32 rounding of inputs and outputs.
sed -i 's/c_float/c_double/g' "$OUT"/xrotor-python/src/*.f90 "$OUT"/xrotor-python/xrotor/*.py
cmake -S "$OUT/xrotor-python" -B "$OUT/xrotor-python/cmake" -DCMAKE_BUILD_TYPE=None \
      -DCMAKE_Fortran_FLAGS="$FFLAGS" > "$OUT/xrotor-python/cmake.log"
cmake --build "$OUT/xrotor-python/cmake" -j 4 > "$OUT/xrotor-python/build.log"
cp "$OUT"/xrotor-python/cmake/libxrotor.so "$OUT/xrotor-python/xrotor/"

# ---------------------------------------------------------------- build record
python3 - "$OUT" "$SRC755" "$SRCXP" "$FFLAGS" <<'EOF'
import json, subprocess, sys
out, s755, sxp, flags = sys.argv[1:]
rev = lambda d: subprocess.run(["git", "-C", d, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
gf = subprocess.run(["gfortran", "--version"], capture_output=True, text=True).stdout.splitlines()[0]
json.dump({"gfortran": gf, "fflags": flags, "xrotor-7.55": rev(s755), "xrotor-python": rev(sxp)},
          open(f"{out}/BUILDINFO.json", "w"), indent=2)
EOF
echo "== done: $OUT"
