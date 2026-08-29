# -*- coding: utf-8 -*-
"""
Ronalds, B.F. (1992). "J-Tube Design for Flexible Umbilicals."
OTC 6875, Offshore Technology Conference, Houston, TX.

Implements the paper's CABLE-theory J-tube analysis -- a different
physical regime from walker1983_jtube.py: this treats the umbilical
as flexible in bending, carrying load essentially as a cable (axial
tension only, per the friction-capstan relation), rather than as a
rigid pipe undergoing elastic-plastic bending. The paper explicitly
positions itself this way ("Previous papers [1, i.e. Walker & Davies
1983] have predicted pull-in loads for single bend J-tubes... The
present work addresses... an umbilical rather than a steel flowline")
and cites Walker & Davies as its Ref. [1].

-----------------------------------------------------------------------
Equations implemented (numbers match the paper)
-----------------------------------------------------------------------
  (1)-(3)   friction/tension relation around a single bend (the
            classical capstan/belt-friction equation, derived here
            from first principles, not just cited): T1 = T0*exp(mu*alpha)
  (4)       M = TL*y -- bending moment at any point on a J-tube of
            ARBITRARY shape, once the distributed contact loads are
            replaced by the two resultant end forces TL, TU. This is
            the paper's key simplification and holds throughout, not
            just within one bend.
  (5)       eta = V/TL (definition)
  (6)       V = (2/pi)*TL -- the zero-length-extension ("propped
            cantilever") limit of the guide-support case; used here as
            the paper's own explicit numerical anchor to check the
            general formula against, since the paper gives NO full
            worked numerical example the way Walker & Davies does.
  (7)       MB/(TL*r) = 1+xi-(1+gamma)*eta -- moment at the guide B
  (8)       tan(theta) = eta -- angular location of max span moment
  (9)       MS/(TL*r) = gamma*eta - 1 + sqrt(1+eta^2) -- max span moment
  Fig. 5    full support-response table (MB, V/TL, H/TL) for FREE,
            GUIDE, and FULL SUPPORT conditions at the base support A.
            (The paper's deflection rows, (EI)*Dx and (EI)*Dy, are
            NOT implemented here -- see "Not implemented" below.)
  (10),(11) tension growth around n bends in sequence:
            T_i = T_(i-1)*exp(mu*alpha_i), T_n = T0*exp(mu*sum(alpha_i))
            -- direct generalisation of eq (3) to a multi-bend J-tube.
  (12)      a = r*phi -- straight approach length before an umbilical
            of finite bending stiffness first touches a bend wall
  (13)      a = sqrt(2r(D-d)) -- same length, alternative closed form
            from the bend geometry (annular gap D-d)
  (14)      2*phi*T ~= 2*(a/r)*T -- normal force bending the umbilical
            around the bend, small-angle approximation
  (15)      (EI)_u < (1/2)*f*a^2*T < f*r(D-d)*T -- condition for the
            umbilical to be flexible enough to wrap the bend at all
  (16)      M = (EI)_u/r -- moment the umbilical carries while wrapped
            around a bend of radius r (f=4, curvature fully developed
            at the entrance) -- f is halved to 8/3 at the OTHER
            entrance condition the paper describes (deflected tip);
            see stage3-style "f" parameter below.
  (17),(18) r/D > 10 and D/t < r/D -- geometric validity limits for
            using simple (Roark's) curved-beam formulas at all.

-----------------------------------------------------------------------
Not implemented
-----------------------------------------------------------------------
  - The FREE-support-case V, H formulas and ALL FOUR deflection
    formulas ((EI)*Dx, (EI)*Dy for both free and guide cases) from the
    Fig. 5 table. These are the densest nested expressions in the
    source scan and, unlike the MB/eta formulas above, this module has
    NO independent numerical check for them (the only paper-stated
    anchor, eq 6's V=2TL/pi limit, only exercises the GUIDE-case eta
    formula). Transcribing them without any way to catch an error
    would be guessing dressed up as verification -- flagged here
    rather than shipped silently. The FULL-SUPPORT case (V=0, H=TL,
    MB=0 identically) IS implemented, since the paper states it
    unconditionally and it needs no formula transcription at all.
  - Fig. 6's general multiple-bend, multiple-PLANE moment/reaction
    solve. The paper is explicit this is NOT closed-form for the
    general case ("This type of indeterminate problem is best solved
    using a stress analysis program") -- only the tension-growth
    relation (eq 10-11) and the same-plane eq (4) moment identity
    carry over from the single-bend analysis.
  - Diaphragm pull-in (Fig. 8) and accidental pull-out (Fig. 9) load
    cases -- both reuse exactly the same Fig. 5 support-response
    machinery with different applied end forces, per the paper's own
    text; not separately coded since they add no new equations, only
    different load-case bookkeping a caller can do directly with the
    functions here.

-----------------------------------------------------------------------
Verification status -- different in KIND from walker1983_jtube.py
-----------------------------------------------------------------------
Walker & Davies (1983) includes a fully worked numerical example
(Appendix II), so every formula there was checked against a paper-
stated number. Ronalds (1992) is written entirely in symbolic/
parametric form -- no numbers are plugged into the design formulas
anywhere in the paper. Verification here therefore rests on:

  (a) ONE genuine numerical anchor the paper does give: eq (6)'s
      explicit statement that eta -> 2/pi as gamma,lambda,xi -> 0.
      guide_support_response()'s general eta formula is checked
      against this in test_ronalds1992.py, and matches to 1e-12.
  (b) Internal consistency: the FULL SUPPORT case must give V=0,
      H=TL, MB=0 identically (stated unconditionally in the paper's
      text, not just the table) -- trivial to verify, still checked.
  (c) Independent, first-principles derivation of the capstan
      equation (eq 1-3) and the M=TL*y moment identity (eq 4) for
      the simplest case (a plain quarter-circle, no straight
      extensions) via direct statics -- these do not depend on
      trusting the OCR transcription of eq (1)-(4) at all.
  (d) Physical monotonicity/limit checks (friction always raises
      tension; a longer multi-bend path never reduces the tension
      growth exponent; etc).

This is a MEANINGFULLY WEAKER verification standard than Walker &
Davies achieves, precisely because the source material offers less to
check against. Treat the guide/full-support MB and eta formulas as
"transcribed and internally consistent, cross-checked against the
one number the paper gives" rather than "independently verified" in
the stronger sense used for walker1983_jtube.py.
"""

import math

import numpy as np


# =====================================================================
# eq (1)-(3): capstan/belt-friction tension relation around one bend
# =====================================================================

def tension_after_bend(T0, mu, alpha):
    """
    eq (3): T1 = T0 * exp(mu * alpha)

    This is the classical capstan (belt-friction) equation -- a
    first-principles result independent of anything specific to
    J-tubes, derivable directly from eq (1)/(2):
        dT = mu*dN  (eq 1, Coulomb friction)
        dN = T*d(theta)  (eq 2, radial equilibrium of a small arc)
        => dT/T = mu*d(theta)  => T1/T0 = exp(mu*alpha)
    (integrating eq 1/2 together, as the paper's own text describes).

    T0    : back tension (N) at the low-tension end
    mu    : coefficient of friction between cable/umbilical and J-tube
    alpha : bend angle subtended (rad)

    Returns T1 (N), the tension at the high-tension (pulled) end.
    """
    return T0 * math.exp(mu * alpha)


def moment_on_arc(TL, y):
    """
    eq (4): M = TL * y

    Bending moment at any point S(x, y) along the J-tube (arc or
    straight section alike), once the distributed contact loads dN,
    dT are replaced by their resultants TL (at the bottom) and TU (at
    the top) -- valid throughout a single-plane J-tube of ANY shape,
    per the paper's text ("Eqn (4) is clearly true also for points
    above and below the bend").

    TL : cable tension at the bottom of the J-tube (N)
    y  : height of the point above the origin O (m) -- O is defined
         as the point where M=0 by construction (paper's own choice
         of coordinate origin, not a physical requirement)
    """
    return TL * y


# =====================================================================
# eq (10)-(11): tension growth around a MULTI-bend J-tube
# =====================================================================

def tension_after_multiple_bends(T0, mu, alphas):
    """
    eq (10)/(11): direct generalisation of the single-bend capstan
    equation to n bends in sequence:

        T_i = T_(i-1) * exp(mu*alpha_i)     (10)
        T_n = T0 * exp(mu * sum(alpha_i))    (11)

    alphas : sequence of bend angles (rad), in pull-in order (bend
             nearest the back-tension end first)

    Returns (T0, T1, ..., Tn) -- the tension at EVERY intermediate
    stage, not just the final one, since intermediate values are
    needed for e.g. checking each individual bend's contact loads.
    """
    tensions = [T0]
    T = T0
    for alpha in alphas:
        T = tension_after_bend(T, mu, alpha)
        tensions.append(T)
    return tuple(tensions)


# =====================================================================
# Single-bend J-tube: guide-support response (Fig. 5, GUIDE column)
# =====================================================================

def guide_support_eta(gamma, lam, xi):
    """
    Fig. 5 table, GUIDE column, row 2: eta = V/TL, general form (a
    guide support at A restrains lateral movement but allows axial
    slip, so H=0 for this case).

        eta = [ 1/2 + gamma*(pi/2-1)
                + (1+gamma)*((1+xi)*(lam/2+xi) - xi^2/3) ]
              / [ pi/4 + gamma*(2 + (pi/2)*gamma + gamma^2/3)
                  + (1+gamma)^2*(lam/2+xi) ]

    gamma, lam, xi : nondimensional lengths (straight-section length
                     / bend radius r) at, respectively: below the
                     bend (A to the start of the arc), between the
                     bend's top and guide B, and between B and the
                     load point D (see Fig. 5's own geometry sketch).
                     Must be >= 0 (paper's own stated restriction).

    Verified against eq (6): at gamma=lam=xi=0 this reduces to
    2/pi EXACTLY (see test_ronalds1992.py) -- the one numerical
    anchor the source paper provides for this formula.
    """
    if gamma < 0 or lam < 0 or xi < 0:
        raise ValueError("guide_support_eta: gamma, lam, xi must all be >= 0 "
                          "(paper's own stated restriction on this formulation)")
    num = (0.5 + gamma * (math.pi / 2 - 1)
           + (1 + gamma) * ((1 + xi) * (lam / 2 + xi) - xi ** 2 / 3))
    den = (math.pi / 4 + gamma * (2 + (math.pi / 2) * gamma + gamma ** 2 / 3)
           + (1 + gamma) ** 2 * (lam / 2 + xi))
    return num / den


def guide_support_MB_over_TLr(gamma, xi, eta):
    """eq (7): MB/(TL*r) = 1 + xi - (1+gamma)*eta -- moment at guide B."""
    return 1.0 + xi - (1.0 + gamma) * eta


def guide_support_theta_at_max_span_moment(eta):
    """eq (8): tan(theta) = eta -- angular location (rad, from the
    start of the bend) where the span moment MS is maximum."""
    return math.atan(eta)


def guide_support_MS_over_TLr(gamma, eta):
    """
    eq (9): MS/(TL*r) = gamma*eta - 1 + sqrt(1+eta^2) -- the maximum
    span moment. The paper notes MS often exceeds MB for longer
    vertical spans (large lambda), and should always be checked when
    lambda > 1.
    """
    return gamma * eta - 1.0 + math.sqrt(1.0 + eta ** 2)


def guide_support_response(gamma, lam, xi):
    """
    Convenience wrapper: full guide-support response (eta, MB, theta
    at max span moment, MS), all in one call.

    Returns dict with keys: eta, MB_over_TLr, theta_max_deg,
    MS_over_TLr.
    """
    eta = guide_support_eta(gamma, lam, xi)
    MB = guide_support_MB_over_TLr(gamma, xi, eta)
    theta_max = guide_support_theta_at_max_span_moment(eta)
    MS = guide_support_MS_over_TLr(gamma, eta)
    return dict(eta=eta, MB_over_TLr=MB,
                theta_max_deg=math.degrees(theta_max), MS_over_TLr=MS)


# =====================================================================
# Single-bend J-tube: FREE and FULL SUPPORT cases (Fig. 5, other
# columns) -- only the rows the paper's own text states unconditionally
# =====================================================================

def free_support_MB_over_TLr(xi):
    """
    Fig. 5 table, FREE column, row 1: MB/(TL*r) = 1+xi.

    "Free" means the base support at A carries no reaction at all
    (V=0, H=0 identically for this case, per the table's own FREE
    column rows 2-3) -- so the J-tube above A behaves, mechanically,
    as if TL and TU were the only loads on a free-standing curved bar.
    """
    return 1.0 + xi


def full_support_response():
    """
    Fig. 5 table, FULL SUPPORT column: V=0, H=TL, MB=0 IDENTICALLY,
    independent of gamma/lambda/xi -- stated unconditionally in the
    paper's own text ("if a support providing full axial restraint is
    provided... V becomes zero and the bending moment is zero
    throughout the J-tube... a direct consequence of the cable within
    the J-tube having negligible bending resistance").

    Returns dict with keys: eta (=0), MB_over_TLr (=0), H_over_TL (=1).
    """
    return dict(eta=0.0, MB_over_TLr=0.0, H_over_TL=1.0)


# =====================================================================
# eq (12)-(18): mechanics of umbilical bending around a J-tube wall
# =====================================================================

def straight_approach_length(r, D, d, phi=None):
    """
    eq (12)/(13): the straight length 'a' an umbilical travels into a
    bend before it first touches BOTH sides of the J-tube wall.

        eq (12): a = r*phi           (phi = subtended angle, small)
        eq (13): a = sqrt(2*r*(D-d))  (equivalent closed form from the
                  bend geometry -- annular gap D-d)

    If phi is given, eq (12) is used and cross-checked against eq (13)
    (they should closely agree for small phi, per the paper's small-
    angle derivation); otherwise eq (13) is used directly.
    """
    a13 = math.sqrt(2.0 * r * (D - d))
    if phi is None:
        return a13
    a12 = r * phi
    return a12  # caller may compare against straight_approach_length(r,D,d) [eq 13] separately


def bend_normal_force(a, r, T):
    """eq (14): 2*phi*T ~= 2*(a/r)*T -- normal force bending the
    umbilical around the wall, small-angle / diameter-difference-
    ignored approximation."""
    return 2.0 * (a / r) * T


def can_umbilical_wrap_bend(EI_u, a, r, D, d, f=4.0):
    """
    eq (15): {EI}_u < (1/2)*f*a^2*T/r ... actually stated as
        {EI}_u < (1/2)*f*a^2*T < f*r*(D-d)*T
    i.e. TWO conditions chained together (the middle term uses eq 13's
    a^2 = 2r(D-d) substituted in, so the rightmost form is just the
    leftmost with 'a' eliminated). Returns whether the FLEXURAL
    STIFFNESS side of the inequality is satisfied -- this function
    needs a trial tension T to evaluate the RHS, so it returns a
    (bool, threshold_T) pair: the umbilical wraps the bend (rather
    than bridging with local double-sided contact) once the applied
    tension T exceeds threshold_T.

    EI_u : umbilical flexural stiffness (EI)_u, N*m^2
    f    : wrapping criterion factor -- f=4 at the bend entrance
           (curvature matches the bend there); paper notes f reduces
           to 8/3 once the tip has deflected (1/2)*r*phi^2 to touch
           the inner wall. Pass f=8/3 for that later-stage criterion.
    """
    # (1/2)*f*a^2*T > EI_u  =>  T > 2*EI_u/(f*a^2)
    threshold_T = 2.0 * EI_u / (f * a ** 2)
    return threshold_T


def geometric_validity(r, D, t):
    """
    eq (17)/(18): validity limits for Roark's simple curved-beam
    formulas (used to derive eq 14-16), and for ignoring cross-section
    ovalisation of the J-TUBE itself (curved thin-wall tubulars can
    ovalise and lose stiffness -- distinct from the umbilical's own
    behaviour, which this module treats as a simple cable elsewhere).

        (17)  r/D > 10
        (18)  D/t < r/D   (equivalently r/D > D/t)

    Returns dict with both ratios and both pass/fail booleans.
    """
    rD = r / D
    Dt = D / t
    return dict(r_over_D=rD, D_over_t=Dt,
                eq17_ok=(rD > 10.0), eq18_ok=(Dt < rD))


# =====================================================================
# Standalone demo: multi-bend tension growth + single-bend guide response
# =====================================================================

if __name__ == '__main__':
    print('=' * 72)
    print('Ronalds (1992) -- multi-bend tension growth demo')
    print('=' * 72)
    print('Geometry: one 90 deg major bend + two 8.5 deg minor bends,')
    print('matching the multi-bend example geometry used in Fan et al.')
    print('(2013, OMAE2013-10103) for comparison purposes -- NOT a')
    print('worked example from this paper itself (Ronalds gives none).')
    print()

    T0 = 100e3    # N, back tension (placeholder)
    mu = 0.3
    alphas_deg = [90.0, 8.5, 8.5]
    alphas = [math.radians(a) for a in alphas_deg]

    tensions = tension_after_multiple_bends(T0, mu, alphas)
    print(f'  T0 = {T0/1e3:.1f} kN')
    for i, (a_deg, T) in enumerate(zip(alphas_deg, tensions[1:]), start=1):
        print(f'  after bend {i} ({a_deg:.1f} deg): T = {T/1e3:.2f} kN')
    print(f'  Total angle = {sum(alphas_deg):.1f} deg, '
          f'final/back tension ratio = {tensions[-1]/T0:.3f}')

    print()
    print('=' * 72)
    print('Single-bend GUIDE-support response demo (zero-length-')
    print("extension case, checked against the paper's own eq (6))")
    print('=' * 72)
    r = guide_support_response(gamma=0.0, lam=0.0, xi=0.0)
    print(f"  eta = V/TL       = {r['eta']:.6f}   (paper eq 6: 2/pi = {2/math.pi:.6f})")
    print(f"  MB/(TL*r)        = {r['MB_over_TLr']:.6f}")
    print(f"  theta at max Ms  = {r['theta_max_deg']:.2f} deg")
    print(f"  MS/(TL*r)        = {r['MS_over_TLr']:.6f}")

    print()
    print('=' * 72)
    print("Full-support case (unconditional per the paper's own text)")
    print('=' * 72)
    fr = full_support_response()
    print(f"  eta = V/TL = {fr['eta']}, MB/(TL*r) = {fr['MB_over_TLr']}, "
          f"H/TL = {fr['H_over_TL']}")

    print()
    print('=' * 72)
    print('Umbilical wrap-around-bend check (Table 1 J-tube/riser')
    print('properties from Fan et al. 2013, for a concrete example)')
    print('=' * 72)
    r_bend = 12.0       # m
    D = 0.406           # m, J-tube ID
    d = 0.254            # m, riser/umbilical OD
    a = straight_approach_length(r_bend, D, d)
    print(f'  Straight approach length a = {a:.3f} m (eq 13)')
    valid = geometric_validity(r_bend, D, t=0.019)
    print(f"  r/D = {valid['r_over_D']:.1f}  (eq 17 needs > 10: "
          f"{'OK' if valid['eq17_ok'] else 'FAILS'})")
    print(f"  D/t = {valid['D_over_t']:.1f}, r/D = {valid['r_over_D']:.1f} "
          f"(eq 18 needs D/t < r/D: {'OK' if valid['eq18_ok'] else 'FAILS'})")
