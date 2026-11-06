"""The transonic lift loss (pyxrotor.transonic) on the P-47 flight tests (validation/naca_p47):
each Clark Y section's lift-divergence Mach number transonic.mlift(t) + shift, with the
drag's Mcrit plain (+0) or calibrated (+0.15). Mean (model - test) efficiency minus the
low-speed bias (tip Mach < 0.9) in tip-Mach bands; writes results.json.

    PYTHONPATH=src python validation/transonic_lift/compare.py
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np

from pyxrotor import Rotor, transonic

HERE = Path(__file__).resolve().parent
VAL = HERE.parent
BANDS = [(0.9, 0.96), (0.96, 1.02), (1.02, 1.08), (1.08, 1.2)]
SETTINGS = {"plain": (0.0, None), "plain, lift loss +0": (0.0, 0.0),
            "calibrated": (0.15, None), "calibrated, lift loss +0": (0.15, 0.0),
            "calibrated, lift loss +0.15": (0.15, 0.15)}


def model(name: str, lift_shift):
    """naca_p47/model.py with its rotors given the lift loss (thickness from its blade form)."""
    spec = importlib.util.spec_from_file_location(name, VAL / "naca_p47" / "model.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    base, thickness = m.sections, {}

    def sections(mcrit_shift=0.0, geo=m.GEO):
        out = base(mcrit_shift, geo)
        r, t = np.array(geo["r_over_R"]), np.array(geo["h_over_b"])
        for sec in out:
            thickness[round(sec.xisect, 9)] = float(np.interp(sec.xisect, r, t))
        return out

    class LiftLossRotor(Rotor):
        def set_sections(self, secs):
            super().set_sections(secs)
            if lift_shift is not None:
                self.set_transonic_lift([transonic.mlift(thickness[round(s.xisect, 9)], lift_shift)
                                         for s in secs])
    m.sections, m.Rotor = sections, LiftLossRotor
    return m


def p47(prop: str, setting: str) -> dict:
    mcrit_shift, lift_shift = SETTINGS[setting]
    m = model(f"tl_p47_{prop}_{setting}", lift_shift)
    res = [(o["Mt"], o["converged"], o["eta"] - float(row["eta"]))
           for row in m.runs(prop) for o in [m.run(row, mcrit_shift=mcrit_shift)]]
    bias = float(np.mean([e for mt, ok, e in res if mt < 0.9 and ok]))
    bands = {}
    for a, b in BANDS:
        e = [e for mt, ok, e in res if a <= mt < b and ok]
        bands[f"{a}-{b}"] = None if not e else round(float(np.mean(e)) - bias, 4)
    return {"bias": round(bias, 4), "band_minus_bias": bands,
            "unconverged": sum(1 for _, ok, _ in res if not ok), "runs": len(res)}


if __name__ == "__main__":
    out = {p: {s: p47(p, s) for s in SETTINGS} for p in ("l177", "l594", "l7d10")}
    (HERE / "results.json").write_text(json.dumps(out, indent=1))
    for p, v in out.items():
        for s, r in v.items():
            print(f"{p:6s} {s:30s} bias {r['bias']:+.3f} {list(r['band_minus_bias'].values())} "
                  f"unconverged {r['unconverged']}/{r['runs']}")
