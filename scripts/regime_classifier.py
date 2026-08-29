# -*- coding: utf-8 -*-
"""
Cross-paper regime classifier: decides whether a given riser/umbilical
belongs to Walker & Davies (1983)'s rigid-pipe, discrete-3-point-
contact regime, or Ronalds (1992)'s flexible-cable, continuous-wrap
regime -- for a given J-tube geometry and (estimated or actual) pull
tension.

Deliberately kept as its own module rather than folded into either
walker1983_jtube.py or ronalds1992_jtube.py: its entire purpose is
comparing across the two papers' domains of applicability, and mixing
it into either paper's own module would blur what that module claims
to implement (see each module's own docstring for how carefully that
provenance is tracked elsewhere in this project).

-----------------------------------------------------------------------
What this answers, and why it matters
-----------------------------------------------------------------------
Walker & Davies and Ronalds are NOT two notations for the same theory
-- they make opposite assumptions about how the pipe/umbilical
contacts the J-tube wall (discrete 3-point contact for a stiff riser
vs continuous wrap-around for a flexible cable), and each paper's
formulas are only valid within its own regime. You cannot use one
paper's numbers to validate the other's method in general.

What CAN be checked, quantitatively: Ronalds' own wrap-vs-bridge
criterion (his eq 15, implemented as
ronalds1992_jtube.can_umbilical_wrap_bend()) gives the tension above
which an object of a given bending stiffness (EI) actually conforms
to the bend wall continuously. Below that threshold it does NOT wrap
-- consistent with Walker's discrete-contact picture instead.

This module computes that threshold for a plain circular tube (given
OD, wall thickness, E) and compares it against an estimated pull
tension, to tell you WHICH paper's assumptions actually apply to your
object -- rather than assuming one or the other.

Verified numerically for Walker's own Appendix II riser (see
test_regime_classifier.py): that riser's wrap threshold is
approximately 12x its own actual Stage I pull load, i.e. nowhere near
Ronalds' regime -- confirming Walker's own worked example is
self-consistent with the discrete-contact assumption it relies on.
"""

import math

import ronalds1992_jtube as rn


def hollow_tube_second_moment(OD, t):
    """
    Second moment of area (m^4) of a thin(ish)-wall circular tube,
    I = (pi/64)*(OD^4 - ID^4), ID = OD - 2*t.

    Standard mechanics-of-materials result, independent of either
    paper -- included here only as a convenience so callers can pass
    plain pipe dimensions (OD, t, E) instead of pre-computing EI
    themselves.
    """
    if t <= 0 or t >= OD / 2.0:
        raise ValueError(f"hollow_tube_second_moment: wall thickness t={t} "
                         f"is not physically possible for OD={OD}")
    ID = OD - 2.0 * t
    return (math.pi / 64.0) * (OD ** 4 - ID ** 4)


def wrap_threshold_tension(OD, t, E, r, D, f=4.0):
    """
    The tension above which an object of outer diameter OD, wall
    thickness t and Young's modulus E would conform continuously to a
    J-tube bend of radius r and inside diameter D (Ronalds' eq 15
    wrap-vs-bridge criterion).

    Returns threshold_T (N). Above this tension: continuous wrap
    (Ronalds' regime). At or below it: discrete/localised contact
    (consistent with Walker's regime, though Walker's own three-point
    model is not itself derived from this criterion -- it is simply
    the OTHER side of the same physical either/or).
    """
    EI = E * hollow_tube_second_moment(OD, t)
    a = rn.straight_approach_length(r, D, OD)
    return rn.can_umbilical_wrap_bend(EI, a, r, D, OD, f=f)


def classify_pullin_regime(OD, t, E, r, D, pullin_tension=None, f=4.0):
    """
    Classifies which paper's assumptions actually apply.

    OD, t, E     : riser/umbilical outer diameter (m), wall thickness
                   (m), Young's modulus (Pa)
    r, D         : J-tube bend radius (m), inside diameter (m)
    pullin_tension : optional estimated/actual pull tension (N). If
                   given, the classification also states the margin
                   (ratio of the wrap threshold to this tension); if
                   omitted, only the threshold itself is returned and
                   no classification label is attempted (a threshold
                   alone doesn't tell you which side of it you're on).
    f            : Ronalds' wrap-criterion factor (default 4.0, "full
                   curvature developed at bend entrance" -- the more
                   conservative of the two values the paper gives; see
                   can_umbilical_wrap_bend()'s docstring for f=8/3).

    Returns dict with keys:
        EI                  : bending stiffness used (N*m^2)
        threshold_T         : wrap-threshold tension (N)
        pullin_tension      : as given (or None)
        margin              : threshold_T / pullin_tension (or None)
        regime              : one of 'rigid (Walker-type)',
                              'flexible (Ronalds-type)',
                              'marginal -- check both', or
                              'unknown (no pull tension given)'
    """
    EI = E * hollow_tube_second_moment(OD, t)
    threshold_T = wrap_threshold_tension(OD, t, E, r, D, f=f)

    if pullin_tension is None:
        return dict(EI=EI, threshold_T=threshold_T, pullin_tension=None,
                    margin=None, regime='unknown (no pull tension given)')

    if pullin_tension <= 0:
        raise ValueError(f"classify_pullin_regime: pullin_tension must be "
                         f"positive, got {pullin_tension}")

    margin = threshold_T / pullin_tension

    # A single hard cutoff at margin=1 would overstate the precision of
    # a criterion that is itself an engineering approximation (Roark's
    # simple beam formulas, small-angle geometry -- see
    # ronalds1992_jtube.py's own "Geometric Limitations" section). A
    # factor-of-2 band either side is treated as "marginal" rather than
    # confidently assigning a regime right at the threshold.
    if margin > 2.0:
        regime = 'rigid (Walker-type): discrete contact, nowhere near wrapping'
    elif margin < 0.5:
        regime = 'flexible (Ronalds-type): well past the wrap threshold'
    else:
        regime = 'marginal -- within 2x of the wrap threshold, check both'

    return dict(EI=EI, threshold_T=threshold_T, pullin_tension=pullin_tension,
                margin=margin, regime=regime)


if __name__ == '__main__':
    print('=' * 72)
    print('Regime classification: Walker & Davies (1983) Appendix II riser')
    print('=' * 72)
    result = classify_pullin_regime(
        OD=0.30, t=0.020, E=200e9,     # Walker Appendix II riser
        r=30.0, D=0.45,                 # Walker Appendix II J-tube
        pullin_tension=167.2e3,          # Walker's own computed P1
    )
    print(f"  EI               = {result['EI']/1e6:.1f} MN*m^2")
    print(f"  Wrap threshold T = {result['threshold_T']/1e6:.2f} MN")
    print(f"  Actual P1        = {result['pullin_tension']/1e3:.1f} kN")
    print(f"  Margin           = {result['margin']:.1f}x")
    print(f"  Regime           = {result['regime']}")

    print()
    print('=' * 72)
    print('Same J-tube, a MUCH more flexible object (placeholder EI ~ a')
    print('typical unbonded flexible pipe/umbilical, order of magnitude')
    print('only -- see flexible_pipe_lib/tier2_stickslip Fig. 6 for real')
    print('measured values on an actual product)')
    print('=' * 72)
    EI_umbilical = 5.0e3   # N*m^2, placeholder -- NOT a real product's stiffness
    threshold = rn.can_umbilical_wrap_bend(
        EI_umbilical,
        rn.straight_approach_length(r=30.0, D=0.45, d=0.30),
        r=30.0, D=0.45, d=0.30)
    print(f"  Placeholder umbilical EI = {EI_umbilical/1e3:.1f} kN*m^2")
    print(f"  Wrap threshold T         = {threshold/1e3:.2f} kN")
    print(f"  (compare: typical umbilical pull tensions are often tens of kN "
          f"and up --")
    print(f"   at that scale this WOULD be past the wrap threshold, "
          f"consistent with")
    print(f"   Ronalds' regime, unlike Walker's steel riser above)")
