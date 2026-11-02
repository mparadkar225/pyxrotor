"""Constants from XROTOR.INC.

PI is deliberately XROTOR's truncated value, not math.pi. Every angle conversion and
radial-grid formula in XROTOR uses it, so using math.pi would put the port about 1e-10
(relative) away from the reference everywhere.
"""

PI = 3.141592654          # PARAMETER (PI=3.141592654)

IX = 100                  # max number of radial stations
IXP = IX + 1
NPARX = 12                # number of stored case parameters
ICASX = 100               # max number of stored cases
IWX = 200                 # work-array dimension
NAX = 20                  # max number of aero sections
NDX = 14                  # parameters per aero section
IQ = IX + 5
JX = (IQ * 3) // 2 + 1
