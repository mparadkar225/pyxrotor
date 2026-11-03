"""How much could slipstream contraction matter for the Report 747 dual propeller?

An estimate, not a model: the rear rotor sees the front's averaged field with radii scaled
by c, axial velocity by 1/c^2 (flux) and swirl by 1/c (circulation). c is the actuator-disk
contraction at the rear's distance, sqrt((1+a)/(1+a(1+z/sqrt(z^2+R^2)))), or fully developed,
sqrt((1+a)/(1+2a)), with a the front's area-weighted mean axial induction factor.

    python validation/naca747/contraction_estimate.py
"""
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "src")]
import numpy as np, model
from pyxrotor.pair import RotorPair, PairConvergenceError

class Contracted(RotorPair):
    c = 1.0
    def _imposed(self, which, other):
        r, u, v = super()._imposed(which, other)
        if which != "rear" or other is None or self.c == 1.0:
            return r, u, v
        uu, _, vv = other.velocity(r / self.c, np.full_like(r, self.spacing))
        return r, uu / self.c**2, self.sign * vv / self.c

def pair(beta, J, c=1.0):
    v, rpm = model.operating_point(J)
    p = Contracted(model.rotor(beta, v), model.rotor(model.REAR[beta], v), model.SPACING, relaxation=0.5)
    p.c = c
    return p, p.solve_rpm(rpm), v

if __name__ == "__main__":
    for beta, J in ((20, 0.4), (25, 0.5), (30, 0.7), (35, 1.2), (40, 1.4), (45, 1.8), (45, 2.2)):
        p0, r0, v = pair(beta, J)
        R = model.R
        rr = np.linspace(model.HUB * R, 0.999 * R, 200)
        a = float(np.sum(p0.front.induced_field().velocity(rr, np.full_like(rr, 1e-4 * R))[0] * rr) / np.sum(rr)) / v
        zr = model.SPACING / R
        for label, c in (("at spacing", np.sqrt((1 + a) / (1 + a * (1 + zr / np.sqrt(zr**2 + 1))))),
                         ("far wake", np.sqrt((1 + a) / (1 + 2 * a)))):
            try:
                _, r1, _ = pair(beta, J, c)
            except PairConvergenceError as e:
                print(beta, J, label, c, "FAILED", e); continue
            print(f"beta {beta} J {J} a={a:.3f} {label:10s} c={c:.4f}: rear power {r1.rear['power']/r0.rear['power']-1:+.2%}, "
                  f"front {r1.front['power']/r0.front['power']-1:+.2%}, rear/front {r0.rear['power']/r0.front['power']:.3f} -> {r1.rear['power']/r1.front['power']:.3f}", flush=True)
