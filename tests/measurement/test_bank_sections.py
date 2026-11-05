"""The transonic propeller comparisons rerun with sections fitted from NASA CR-176883's data
banks (validation/cr176883_sections) instead of pyxrotor.sections' thickness estimates:

  * the bank sections alone (their own fitted Mcrit) still lose far too much at transonic
    tips: on the Clark Y propellers (P-47 x3, YP-39) 40-125% of 7.55's excess loss at tip
    Mach 1.02-1.08 (the bank's Mcrit is up to 0.07 above the estimate at 5-8% but not at
    the 4% tips); on the 16-series ones (YP-47M, SR-2, L51L28) *more* than 7.55, because the
    16-series bank's drag-rise Mach number is 0.03-0.07 below the estimate at most design CLs
    (and SR-2's low design CLs use the bank's CLD 0.2);
  * they need a shift about as large as the estimates: Clark Y + 0.10-0.15 (P-47) and
    + 0.15-0.20 (YP-39); 16-series + 0.15-0.20 (YP-47M 836), + 0.20-0.25 (L51L28, SR-2).
    So the calibration is not a property of the section estimate it was made on;
  * the low-speed level barely changes on the Clark Y propellers (P-47 bias within 0.02 of
    the estimates'). The 16-series bank's larger drag due to lift raises the static power of
    the L51L28 blades by 6-14%, which narrows the static under-prediction on the cambered
    blades (CP at design CL 0.3, 4 and 8 deg: 0.94/0.73 -> 1.04/0.78 of the test's).
"""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

pytestmark = pytest.mark.measurement
VAL = Path(__file__).resolve().parents[2] / "validation"
spec = importlib.util.spec_from_file_location("bank_sections_compare",
                                              VAL / "cr176883_sections" / "compare.py")
cmp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cmp)


def worst(d):
    return max(abs(v) for v in d.values() if v is not None)


@pytest.mark.parametrize("prop,shift", [("l177", "bank+0.15"), ("l594", "bank+0.10"),
                                        ("l7d10", "bank+0.15")])
def test_p47(prop, shift):
    bank0, best, est = (cmp.p47(prop, s) for s in ("bank+0.00", shift, "est+cal"))
    assert worst(best["band_minus_bias"]) < 0.07 and best["unconverged"] <= 2
    assert min(v for v in bank0["band_minus_bias"].values() if v is not None) < -0.2
    assert abs(best["bias"] - est["bias"]) < 0.02


def test_yp39():
    assert worst(cmp.yp39("bank+0.15")["band"]) < 0.04
    assert cmp.yp39("bank+0.00")["band"]["1.035-1.07"] < -0.4


def test_yp47m():
    r0, r15, r20 = (cmp.yp47m(s) for s in ("bank+0.00", "bank+0.15", "bank+0.20"))
    e0 = cmp.yp47m("est+0")
    keys = ("M0 0.3 J 1.0", "M0 0.5 J 1.75")
    assert all(r0[k] < e0[k] for k in keys)          # bank alone loses more than 7.55
    assert all(abs(r15[k]) < 0.04 and abs(r20[k]) < 0.02 for k in keys)


def test_sr2():
    r = {s: cmp.sr2(s)["minus_bias"] for s in ("est+0", "bank+0.00", "bank+0.15", "bank+0.25")}
    assert r["bank+0.00"]["0.85"] < r["est+0"]["0.85"] < -0.5
    assert r["bank+0.15"]["0.85"] < -0.2
    assert worst(r["bank+0.25"]) < 0.05


def test_l51l28():
    b20, est = cmp.l51l28("bank+0.20"), cmp.l51l28("est+cal")
    for k in ("cld 0.3 beta 4.0", "cld 0.3 beta 8.0", "cld 0.3 beta 12.0", "cld 0.5 beta 8.0",
              "cld 0.5 beta 12.0"):
        assert abs(b20[k]["drift_Mt1.0"] - 1.0) < 0.1, k
    for k in ("cld 0.3 beta 4.0", "cld 0.3 beta 8.0", "cld 0.5 beta 4.0", "cld 0.5 beta 8.0"):
        ratio = b20[k]["CP_ratio_Mt0.4"] / est[k]["CP_ratio_Mt0.4"]
        assert 1.05 < ratio < 1.17, (k, ratio)
    assert np.isclose(b20["cld 0.3 beta 4.0"]["CP_ratio_Mt0.4"], 1.04, atol=0.01)
