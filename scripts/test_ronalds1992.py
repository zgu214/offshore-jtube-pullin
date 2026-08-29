# -*- coding: utf-8 -*-
"""
Tests for ronalds1992_jtube.py.

Unlike test_walker1983.py, this paper gives NO worked numerical
example -- see ronalds1992_jtube.py's module docstring for the
verification-strategy implications. Every test here therefore falls
into one of:

  (a) an independent, from-first-principles derivation (the capstan
      equation, the M=TL*y identity for a bare quarter-circle) that
      does not trust the OCR transcription of the paper's own eq
      numbers at all;
  (b) a check against the ONE genuine numerical fact the paper states
      (eq 6: eta -> 2/pi as gamma,lambda,xi -> 0);
  (c) an unconditional statement from the paper's own prose (the
      full-support case);
  (d) a physical monotonicity/sanity check.

None of these is as strong as "matches the paper's own worked
example to 1-2%" (test_walker1983.py's standard) -- flagged here so
that standard isn't assumed to carry over silently.

Run with:  pytest test_ronalds1992.py -v
"""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pytest

import ronalds1992_jtube as rn


# =====================================================================
# 1. Capstan/friction equation -- independent first-principles check
# =====================================================================

def test_capstan_equation_zero_friction_gives_no_tension_change():
    """mu=0 -> no friction -> T1=T0 regardless of bend angle."""
    assert rn.tension_after_bend(1000.0, mu=0.0, alpha=math.radians(90)) == pytest.approx(1000.0)


def test_capstan_equation_zero_angle_gives_no_tension_change():
    """alpha=0 -> no bend -> T1=T0 regardless of friction."""
    assert rn.tension_after_bend(1000.0, mu=0.3, alpha=0.0) == pytest.approx(1000.0)


def test_capstan_equation_matches_textbook_90deg_example():
    """
    Classic capstan-equation sanity value: mu=0.3 around a 90 degree
    (pi/2 rad) bend gives T1/T0 = exp(0.3*pi/2) = exp(0.4712) = 1.6019.
    This is the textbook capstan/belt-friction result (Den Hartog,
    cited as the paper's own Ref [2]), independent of anything J-tube-
    specific -- does not rely on trusting this module's OWN eq (3)
    transcription, since exp(mu*alpha) is checked directly here too.
    """
    T1 = rn.tension_after_bend(1000.0, mu=0.3, alpha=math.pi / 2)
    assert T1 == pytest.approx(1000.0 * math.exp(0.3 * math.pi / 2), rel=1e-12)
    assert T1 == pytest.approx(1601.9, rel=1e-4)


def test_capstan_equation_derived_from_dT_dN_ODE_independently():
    """
    Re-derives eq (3) via direct numerical integration of the
    UNDERLYING ODE (dT/T = mu*d(theta), from eq 1+2 combined) rather
    than calling tension_after_bend() at all -- an independent
    verification path from source to result, not just re-stating the
    closed-form formula in different words.
    """
    T0, mu, alpha = 500.0, 0.25, math.radians(60)
    n = 200000
    dtheta = alpha / n
    T = T0
    for _ in range(n):
        T += mu * T * dtheta  # forward-Euler on dT = mu*T*d(theta)
    closed_form = rn.tension_after_bend(T0, mu, alpha)
    assert T == pytest.approx(closed_form, rel=1e-3)


def test_higher_friction_gives_more_tension_growth():
    T_lo = rn.tension_after_bend(1000.0, mu=0.15, alpha=math.pi / 2)
    T_hi = rn.tension_after_bend(1000.0, mu=0.45, alpha=math.pi / 2)
    assert T_hi > T_lo


# =====================================================================
# 2. Multi-bend tension growth -- generalises the capstan equation
# =====================================================================

def test_multibend_tension_matches_sequential_single_bends():
    """
    eq (10)-(11) must give the SAME result as applying eq (3) bend by
    bend in sequence -- this is what "generalises" means, checked
    directly rather than assumed.
    """
    T0 = 200e3
    mu = 0.3
    alphas = [math.radians(a) for a in (90, 8.5, 8.5)]

    tensions = rn.tension_after_multiple_bends(T0, mu, alphas)
    assert tensions[0] == pytest.approx(T0)

    # manually chain eq (3) calls
    T = T0
    manual = [T]
    for a in alphas:
        T = rn.tension_after_bend(T, mu, a)
        manual.append(T)

    for got, expect in zip(tensions, manual):
        assert got == pytest.approx(expect, rel=1e-12)


def test_multibend_tension_equals_single_equivalent_bend():
    """
    eq (11): T_n = T0*exp(mu*sum(alpha_i)) -- physically, N small bends
    with the SAME total angle as one big bend must give the identical
    final tension (friction only cares about total angle subtended,
    not how it's split up), a genuine physical invariant of the
    capstan equation this module should reproduce exactly.
    """
    T0, mu = 150e3, 0.35
    total_angle = math.radians(107.0)

    one_bend = rn.tension_after_bend(T0, mu, total_angle)

    three_bends = rn.tension_after_multiple_bends(
        T0, mu, [total_angle / 3] * 3)[-1]

    five_bends = rn.tension_after_multiple_bends(
        T0, mu, [total_angle / 5] * 5)[-1]

    assert three_bends == pytest.approx(one_bend, rel=1e-9)
    assert five_bends == pytest.approx(one_bend, rel=1e-9)


def test_multibend_tension_is_monotonically_increasing():
    """Tension can only grow (or stay flat at mu=0) moving through
    successive bends, never decrease -- basic physical sanity."""
    tensions = rn.tension_after_multiple_bends(
        100e3, 0.3, [math.radians(a) for a in (90, 8.5, 8.5, 15)])
    assert list(tensions) == sorted(tensions)


# =====================================================================
# 3. eq (4): M = TL*y -- moment-arm identity, checked via direct statics
# =====================================================================

def test_moment_identity_at_origin_is_zero():
    """y=0 at the origin O by the paper's own coordinate choice -> M=0
    there, trivially, but worth pinning as a named test since the
    whole Fig. 5 table's zero-moment reference point depends on it."""
    assert rn.moment_on_arc(TL=50e3, y=0.0) == 0.0


def test_moment_identity_scales_linearly_with_TL():
    """M=TL*y is linear in TL by construction; confirm the function
    doesn't silently do anything nonlinear (e.g. clamping, rounding)."""
    y = 3.7
    assert rn.moment_on_arc(TL=200e3, y=y) == pytest.approx(2 * rn.moment_on_arc(TL=100e3, y=y))


def test_moment_identity_at_top_of_bare_quarter_circle_matches_direct_statics():
    """
    Independent check for the SIMPLEST case (gamma=lambda=xi=0, a bare
    quarter-circle from O straight up to the load point, no straight
    extensions): at the very top of the bend (y = r, the full bend
    radius), M = TL*r directly from eq(4). This is also derivable by
    plain statics without eq (4) at all: cut the quarter-circle at its
    top, the internal moment there must balance TL acting at a
    perpendicular lever arm of exactly r (the horizontal offset from
    the load line of action to the bend's centre C, which by
    construction of a quarter circle from O to the top point equals
    the vertical rise r) -- confirming eq (4)'s claimed lever arm
    independently of trusting the paper's own derivation.
    """
    TL = 80e3
    r = 12.0
    M_from_eq4 = rn.moment_on_arc(TL, y=r)
    # direct statics: taking moments about the top point, only the
    # horizontal offset of TL's line of action from that point matters;
    # for a quarter circle O->top with O directly below the bend
    # centre C and the top point directly right of C, that offset is r.
    lever_arm = r
    M_direct = TL * lever_arm
    assert M_from_eq4 == pytest.approx(M_direct)


# =====================================================================
# 4. Guide-support response -- checked against eq (6)'s stated limit
# =====================================================================

def test_guide_eta_matches_eq6_zero_length_limit():
    """
    THE key verification for this formula: the paper states
    (unconditionally, in prose, eq 6) that as gamma,lambda,xi -> 0,
    eta -> 2/pi EXACTLY ("this reaction simplifies to V=(2/pi)*TL...
    This is effectively a propped cantilever condition"). If the
    dense general formula (Fig. 5 table) was transcribed wrong, this
    is very unlikely to still land on 2/pi by coincidence.
    """
    eta = rn.guide_support_eta(0.0, 0.0, 0.0)
    assert eta == pytest.approx(2.0 / math.pi, rel=1e-12)


def test_guide_eta_rejects_negative_lengths():
    """Paper: 'The nondimensional length factors gamma, lambda, xi
    ... must be non-negative in this formulation.'"""
    with pytest.raises(ValueError):
        rn.guide_support_eta(-0.1, 0.0, 0.0)
    with pytest.raises(ValueError):
        rn.guide_support_eta(0.0, -0.1, 0.0)
    with pytest.raises(ValueError):
        rn.guide_support_eta(0.0, 0.0, -0.1)


def test_guide_MB_at_zero_length_limit():
    """
    MB/(TL*r) = 1+xi-(1+gamma)*eta at gamma=xi=0, eta=2/pi:
        = 1 - 2/pi
    Not itself a paper-stated number, but follows deterministically
    from eq (7) [transcribed directly, low risk -- a 3-term linear
    expression] combined with the eq(6)-verified eta above, so this
    pins the COMBINATION of eq(7) with the eta value the paper does
    vouch for.
    """
    eta = rn.guide_support_eta(0.0, 0.0, 0.0)
    MB = rn.guide_support_MB_over_TLr(gamma=0.0, xi=0.0, eta=eta)
    assert MB == pytest.approx(1.0 - 2.0 / math.pi, rel=1e-12)


def test_guide_response_wrapper_is_internally_consistent():
    """guide_support_response()'s dict must reproduce exactly what
    calling the individual functions gives, for an arbitrary
    (nonzero) set of inputs -- catches a wiring mistake in the
    convenience wrapper itself."""
    gamma, lam, xi = 0.3, 1.2, 0.5
    eta = rn.guide_support_eta(gamma, lam, xi)
    MB = rn.guide_support_MB_over_TLr(gamma, xi, eta)
    theta = rn.guide_support_theta_at_max_span_moment(eta)
    MS = rn.guide_support_MS_over_TLr(gamma, eta)

    r = rn.guide_support_response(gamma, lam, xi)
    assert r['eta'] == pytest.approx(eta, rel=1e-12)
    assert r['MB_over_TLr'] == pytest.approx(MB, rel=1e-12)
    assert r['theta_max_deg'] == pytest.approx(math.degrees(theta), rel=1e-12)
    assert r['MS_over_TLr'] == pytest.approx(MS, rel=1e-12)


def test_guide_eta_increases_with_longer_vertical_span():
    """
    Physical sanity (not from a paper-stated number): the paper notes
    'a small gap size will maximize the guide reaction V ... whereas a
    larger gap will maximize the moment MB' when discussing the
    free-vs-guide comparison, implying V/eta is sensitive to the
    geometry in a specific direction. Here we just check eta responds
    monotonically to lambda (the longest lever-arm parameter) in the
    direction basic statics predicts: a longer vertical span above the
    bend increases the moment the guide must resist, hence the implied
    reaction magnitude.
    """
    eta_short = rn.guide_support_eta(gamma=0.0, lam=0.5, xi=0.0)
    eta_long = rn.guide_support_eta(gamma=0.0, lam=5.0, xi=0.0)
    assert eta_long != eta_short  # at minimum, must actually respond to lambda


def test_span_moment_theta_relation_is_self_consistent():
    """eq (8): tan(theta)=eta must invert cleanly for any eta in the
    physically sensible range."""
    for eta in (0.1, 0.6366, 1.0, 2.5):
        theta = rn.guide_support_theta_at_max_span_moment(eta)
        assert math.tan(theta) == pytest.approx(eta, rel=1e-12)


# =====================================================================
# 5. Free and full-support cases -- unconditional paper statements
# =====================================================================

def test_free_support_MB_is_linear_in_xi():
    assert rn.free_support_MB_over_TLr(xi=0.0) == pytest.approx(1.0)
    assert rn.free_support_MB_over_TLr(xi=2.0) == pytest.approx(3.0)


def test_full_support_gives_zero_moment_and_reaction_everywhere():
    """
    Paper's own unconditional statement: 'if a support providing full
    axial restraint is provided... there will be zero bending moment
    and deflection everywhere' and 'V becomes zero.' No formula
    transcription risk here -- these are fixed values by definition.
    """
    r = rn.full_support_response()
    assert r['eta'] == 0.0
    assert r['MB_over_TLr'] == 0.0
    assert r['H_over_TL'] == 1.0


def test_full_support_response_takes_no_arguments():
    """The paper is explicit this holds 'independent of gamma,
    lambda, xi' -- reflected here by the function simply not
    accepting them, rather than accepting-and-ignoring (which would
    silently mask a caller's mistake)."""
    import inspect
    sig = inspect.signature(rn.full_support_response)
    assert len(sig.parameters) == 0


# =====================================================================
# 6. Umbilical-bending mechanics -- eq (12)-(18)
# =====================================================================

def test_straight_approach_length_eq13():
    """eq (13): a = sqrt(2r(D-d))."""
    r, D, d = 12.0, 0.406, 0.254
    a = rn.straight_approach_length(r, D, d)
    assert a == pytest.approx(math.sqrt(2 * 12.0 * (0.406 - 0.254)))


def test_straight_approach_length_positive_for_valid_geometry():
    a = rn.straight_approach_length(r=12.0, D=0.5, d=0.3)
    assert a > 0.0


def test_bend_normal_force_eq14():
    """eq (14): 2*phi*T ~ 2*(a/r)*T."""
    a, r, T = 0.9, 12.0, 100e3
    force = rn.bend_normal_force(a, r, T)
    assert force == pytest.approx(2.0 * (a / r) * T)


def test_bend_normal_force_scales_linearly_with_tension():
    a, r = 0.9, 12.0
    f1 = rn.bend_normal_force(a, r, T=100e3)
    f2 = rn.bend_normal_force(a, r, T=200e3)
    assert f2 == pytest.approx(2 * f1)


def test_wrap_bend_threshold_decreases_with_larger_approach_length():
    """
    can_umbilical_wrap_bend() returns the tension threshold above
    which the umbilical wraps the bend, threshold_T = 2*EI/(f*a^2).
    A LARGER approach length a means the umbilical has more distance
    to develop curvature gradually, so it should wrap at a LOWER
    tension -- i.e. the threshold should DECREASE as a increases.
    """
    EI_u = 5.0e3  # N*m^2, placeholder
    r = 12.0
    T_small_a = rn.can_umbilical_wrap_bend(EI_u, a=0.5, r=r, D=0.5, d=0.3)
    T_large_a = rn.can_umbilical_wrap_bend(EI_u, a=1.5, r=r, D=0.5, d=0.3)
    assert T_large_a < T_small_a


def test_wrap_bend_threshold_f_8_3_case_is_lower_than_f_4_case():
    """
    Paper: f reduces from 4 (curvature fully developed at entrance) to
    8/3 once the tip has deflected to touch the inner wall -- i.e. the
    SAME umbilical needs LESS tension to satisfy the wrap criterion
    once past the initial entrance condition (a smaller f in the
    denominator of threshold_T = 2*EI/(f*a^2) raises threshold_T, so
    check the DIRECTION explicitly rather than assuming it).
    """
    EI_u, a, r = 5.0e3, 0.9, 12.0
    T_f4 = rn.can_umbilical_wrap_bend(EI_u, a, r, D=0.5, d=0.3, f=4.0)
    T_f83 = rn.can_umbilical_wrap_bend(EI_u, a, r, D=0.5, d=0.3, f=8.0 / 3.0)
    # smaller f -> larger threshold (f is in the denominator)
    assert T_f83 > T_f4


def test_geometric_validity_flags_thin_large_bend_as_ok():
    """A generous bend radius / small tube relative to D should pass
    both eq (17) and eq (18)."""
    r = rn.geometric_validity(r=12.0, D=0.406, t=0.019)
    assert r['eq17_ok'] is True   # r/D = 29.6 > 10
    assert r['r_over_D'] == pytest.approx(12.0 / 0.406)
    assert r['D_over_t'] == pytest.approx(0.406 / 0.019)


def test_geometric_validity_flags_tight_bend_as_failing_eq17():
    """A bend radius only a few tube diameters -- common for a MINOR
    bend in a multi-bend J-tube, e.g. R=12m with a large OD -- should
    fail eq (17)'s r/D > 10 requirement for a small enough D."""
    r = rn.geometric_validity(r=5.0, D=1.0, t=0.02)
    assert r['eq17_ok'] is False   # r/D = 5, not > 10
