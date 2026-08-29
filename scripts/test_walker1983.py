# -*- coding: utf-8 -*-
"""
Tests for walker1983_jtube.py.

Two kinds of test here:

1. Structural/limit checks on the building blocks (moment-curvature,
   alpha/beta), independent of any specific worked example -- these
   catch a broken formula even if the Appendix II numbers happened to
   still line up by coincidence.

2. Reproduction of the Appendix II worked example, at the tolerances
   documented in walker1983_jtube.py's module docstring (1-2%, from
   the paper reading Fig. 12 by eye rather than evaluating its own
   equations). These tolerances are set FROM that documented
   discrepancy, not loosened after the fact to force a pass.

Run with:  pytest test_walker1983.py -v
"""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pytest

import walker1983_jtube as jt


# =====================================================================
# Appendix II reference case
# =====================================================================

SIGMA_O = 400e6   # Pa
E = 200e9         # Pa
R = 30.0          # m  (R/d = 100)
D = 0.45          # m  (D/d = 1.5)
d = 0.30          # m
t = 0.020         # m
MU = 0.30


@pytest.fixture(scope="module")
def result():
    return jt.full_pullin_analysis(R, D, d, t, SIGMA_O, E, MU)


# =====================================================================
# 1. Building blocks -- structural checks, no worked example needed
# =====================================================================

def test_Mo_is_pi_over_4_times_Mp():
    """eq (4): Mo = (pi/4) Mp is an exact algebraic identity."""
    Mp, Mo = jt.Mp_Mo(dm=0.28, t=0.02, sigma_o=400e6)
    assert Mo == pytest.approx(Mp * math.pi / 4, rel=1e-12)


def test_moment_ratio_quadrant_endpoints():
    """
    eq (3) at its two defining limits:
      omega=0    -> M/(dm^2 t sigma_o) = 1        (= Mo/(dm^2 t sigma_o), first yield)
      omega=pi/2 -> M/(dm^2 t sigma_o) = pi/4 ... wait, defined as Mp -- see below.
    Actually eq(3)'s own normalization is by (dm^2*t*sigma_o) = Mp, so:
      omega=0    -> sin(0)+sec(0)*(pi/2-0) = 0 + 1*(pi/2) = pi/2 ... NOT 1.
    Re-derive directly from the formula instead of assuming the answer.
    """
    at_0 = jt.moment_ratio_quadrant(np.array([0.0]))[0]
    at_max = jt.moment_ratio_quadrant(np.array([math.pi / 2]))[0]
    # omega=0: sin(0) + sec(0)*(pi/2 - 0) = 0 + 1*(pi/2) = pi/2
    assert at_0 == pytest.approx(math.pi / 2, rel=1e-9)
    # omega=pi/2 is the removable singularity of sec(omega)*(pi/2-omega);
    # the function's own limit is exactly 1 (fully plastic, M=Mp).
    assert at_max == pytest.approx(1.0, rel=1e-9)


def test_alpha_beta_equal_one_at_k_equals_1():
    """
    eq (1): Delta=alpha*Delta_e, theta=beta*theta_e. At k=1 (M=Mo, the
    first-yield boundary) there is NO plastic straining yet, so the
    cantilever response must be exactly elastic: alpha=beta=1.
    """
    assert jt.alpha_factor(1.0) == pytest.approx(1.0, abs=1e-12)
    assert jt.beta_factor(1.0) == pytest.approx(1.0, abs=1e-12)


def test_alpha_beta_are_one_below_k_equals_1():
    """Below first yield the response is elastic by definition (eq 1),
    not by extrapolating the eq (7)/(8) polynomial fits."""
    for k in (0.1, 0.5, 0.9, 0.999):
        assert jt.alpha_factor(k) == 1.0
        assert jt.beta_factor(k) == 1.0


def test_alpha_beta_increase_with_k_above_1():
    """
    Fig. 12: both curves rise monotonically (steeply) as k increases
    past 1 -- more plastic straining means more deflection/rotation
    per unit root moment relative to the elastic prediction.
    """
    ks = np.array([1.0, 1.05, 1.1, 1.2, 1.3])
    a = jt.alpha_factor(ks)
    b = jt.beta_factor(ks)
    assert np.all(np.diff(a) > 0)
    assert np.all(np.diff(b) > 0)


def test_strain_ratio_elastic_plastic_branches_agree_at_boundary():
    """eq (5) and eq (6) must meet exactly at M/Mo=1, eps/eps_o=1."""
    from_plastic = jt.strain_ratio_from_moment(1.0, branch='plastic')
    from_elastic = jt.strain_ratio_from_moment(1.0, branch='elastic')
    assert from_plastic == pytest.approx(1.0, abs=1e-9)
    assert from_elastic == pytest.approx(1.0, abs=1e-9)


def test_strain_ratio_elastic_branch_is_linear_in_MMo():
    """
    eq (6) is (4/pi)*(M/Mp) in the PAPER's own variable; converted to
    this function's M/Mo argument it collapses to the trivial identity
    eps/eps_o = M/Mo (strain proportional to moment below first
    yield). Catches the specific bug where the (4/pi) factor was left
    in after converting M/Mp -> M/Mo, which would break the M/Mo=1
    boundary-agreement test above only by coincidence of testing
    exactly at the one point where both forms agree (MMo=1).
    """
    for MMo in (0.2, 0.5, 0.8, 0.99):
        assert jt.strain_ratio_from_moment(MMo, branch='elastic') == pytest.approx(MMo, rel=1e-12)


def test_strain_ratio_inversion_round_trips():
    """MMo_from_strain_ratio() must invert strain_ratio_from_moment()."""
    for k_true in (1.05, 1.1, 1.2351, 1.3):
        eps_ratio = jt.strain_ratio_from_moment(k_true)
        k_back = jt.MMo_from_strain_ratio(eps_ratio)
        assert k_back == pytest.approx(k_true, rel=1e-6)


def test_k_from_alpha_Mbar_Mp_round_trips():
    """Stage III's inverse lookup must invert alpha(k)*(k*pi/4)."""
    for k_true in (1.0, 1.05, 1.2, 1.3):
        target = jt.alpha_factor(k_true) * (k_true * math.pi / 4.0)
        k_back = jt.k_from_alpha_Mbar_Mp(target)
        assert k_back == pytest.approx(k_true, rel=1e-6)


def test_k_from_alpha_Mbar_Mp_below_elastic_threshold():
    """
    Below pi/4 (=alpha*Mbar/Mp at k=1), alpha=1 identically, so
    k = target*(4/pi) directly -- no root-find needed, and the direct
    formula must agree with what a root-find on the (trivial) elastic
    branch would give.
    """
    target = 0.5   # < pi/4 = 0.7854
    k = jt.k_from_alpha_Mbar_Mp(target)
    assert k == pytest.approx(target * 4.0 / math.pi, rel=1e-12)
    assert k < 1.0


def test_pullin_load_ratio_prefactor_multiplies_whole_bracket():
    """
    eq (19): confirms (Mbar/pi/Mp) multiplies BOTH the plasticity term
    (d/R - 2 sigma_o/E) and the friction term 2*mu*(d/L+K) -- easy to
    mis-parse as only multiplying the first term. Construct a case
    where getting this wrong changes the answer by a large factor, so
    the test would fail loudly if the prefactor were misapplied.
    """
    Mbar_Mp = 0.5
    Rd = 100.0
    sigma_o_E = 0.002
    mu = 0.3
    Ld = 30.0
    Kd = 0.05
    got = jt.pullin_load_ratio(Mbar_Mp, Rd, sigma_o_E, mu, Ld, Kd)
    bracket = (1.0 / Rd) - 2 * sigma_o_E + 2 * mu * (1.0 / Ld + Kd)
    expected = (Mbar_Mp / math.pi) * bracket
    assert got == pytest.approx(expected, rel=1e-12)
    # and NOT the "prefactor only on the plasticity term" misreading:
    wrong = (Mbar_Mp / math.pi) * (1.0 / Rd - 2 * sigma_o_E) + 2 * mu * (1.0 / Ld + Kd)
    assert got != pytest.approx(wrong, rel=1e-3)


# =====================================================================
# 2. Appendix II worked example -- end-to-end reproduction
# =====================================================================

def test_Po_matches_paper_exactly(result):
    """
    Po = pi*d*t*sigma_o (outer diameter). This is the one place the
    module docstring flags a genuine inconsistency in the source paper
    (Fig. 7 caption vs. the worked example) -- see full_pullin_analysis().
    The worked example's own number is reproduced to <0.01%.
    """
    assert result['Po'] / 1e3 == pytest.approx(7540.0, rel=1e-3)


def test_stage1_Mbar_Mp(result):
    assert result['stage1']['Mbar_Mp'] == pytest.approx(0.97, rel=5e-3)


def test_stage1_L_over_d(result):
    assert result['stage1']['L_over_d'] == pytest.approx(32.7, rel=0.02)


def test_stage1_l1_over_d(result):
    assert result['stage1']['l1_over_d'] == pytest.approx(12.6, rel=0.02)


def test_stage1_tau(result):
    assert result['stage1']['tau_deg'] == pytest.approx(7.1, rel=0.03)


def test_stage1_P1_over_Po(result):
    assert result['stage1']['P1_over_Po'] == pytest.approx(2.2e-2, rel=0.02)


def test_stage1_P1_kN(result):
    assert result['P1'] / 1e3 == pytest.approx(166.0, rel=0.02)


def test_direct_entry_case_raises_rather_than_silently_matching(result):
    """
    Fig. 13's direct-entry geometry is NOT implemented (see
    stage1_l1_over_d_direct_entry()'s docstring: two straightforward
    readings of the paper's one-paragraph description were tried and
    both missed the paper's own stated result by more than this
    module's usual 1-2% tolerance). It must raise, not silently return
    the standard-case answer unchanged -- a caller asking for
    direct_entry=True and getting the SAME P1 as the lead-in case would
    be a silent, undetectable wrong answer, which is worse than an
    explicit NotImplementedError.
    """
    with pytest.raises(NotImplementedError):
        jt.full_pullin_analysis(R, D, d, t, SIGMA_O, E, MU, direct_entry=True)


def test_stage2_C1_C2(result):
    assert result['stage2']['C1'] == pytest.approx(125.0, rel=1e-3)
    assert result['stage2']['C2'] == pytest.approx(115.4, rel=2e-3)


def test_stage2_l2_over_d(result):
    assert result['stage2']['l2_over_d'] == pytest.approx(34.5, rel=0.01)


def test_stage2_P2_over_Po(result):
    assert result['stage2']['P2_over_Po'] == pytest.approx(1.27e-2, rel=0.03)


def test_stage2_C1_C2_no_plastic_straining_variant_is_smaller():
    """
    The paper gives an alternative C1=R/2d, C2=2R/3d for 'the
    situation where no plastic straining has occurred'. For a bend
    radius that DOES cause plastic straining (as here), that
    simplified pair should differ from (and be smaller than) the
    general plastic-straining pair -- otherwise the two formulas
    would be indistinguishable and the paper wouldn't bother giving
    both.
    """
    Rd, sigma_o_E = 100.0, 0.002
    C1_gen, C2_gen = jt.stage2_C1_C2(Rd, sigma_o_E, no_plastic_straining=False)
    C1_np, C2_np = jt.stage2_C1_C2(Rd, sigma_o_E, no_plastic_straining=True)
    assert C1_np == pytest.approx(Rd / 2.0)
    assert C2_np == pytest.approx(2.0 * Rd / 3.0)
    assert C1_np < C1_gen
    assert C2_np < C2_gen


def test_stage3_Rp_over_d(result):
    assert result['stage3']['Rp_over_d'] == pytest.approx(166.7, rel=2e-3)


def test_stage3_peak_near_papers_middle_trial(result):
    """
    The paper explores l3/d in {17.5, 20, 22.5} by hand and finds
    l3/d=20 gives the best (highest) P'/Po of the three (0.0092).
    This module scans a fine grid and finds the true peak, which
    should land close to that same middle trial (not at 17.5 or 22.5,
    confirming the paper's hand search wasn't just lucky/coarse-grained
    in a way that missed a materially different optimum).
    """
    s3 = result['stage3']
    assert 15.0 < s3['l3_over_d'] < 25.0
    assert s3['P_prime_over_Po'] == pytest.approx(0.0092, rel=0.05)


def test_stage3_P_prime_is_the_scan_maximum(result):
    """The reported peak must actually be the max of the returned scan,
    not just some interior point."""
    scan = result['stage3']['scan']
    assert result['stage3']['P_prime_over_Po'] == pytest.approx(
        scan[:, 4].max(), rel=1e-9)


def test_stage3_P3_over_Po(result):
    assert result['stage3']['P3_over_Po'] == pytest.approx(2.19e-2, rel=0.03)


def test_P3_greater_than_P2_less_than_P1(result):
    """
    Physical ordering from the paper's Fig. 3 test record: Stage II is
    a LOCAL MINIMUM, Stage III rises back to about the same order as
    Stage I ('the approximate equality of the Stage I and Stage III
    pull-in loads has been observed experimentally for a wide range').
    """
    P1, P2, P3 = result['P1'], result['P2'], result['P3']
    assert P2 < P1
    assert P2 < P3
    assert 0.8 < P3 / P1 < 1.2   # "of the same order", per the paper


# =====================================================================
# 3. Sanity on parameter sensitivity (no reference numbers -- just
#    checking the model responds the way basic mechanics demands)
# =====================================================================

def test_higher_friction_gives_higher_pullin_load():
    r_lo = jt.full_pullin_analysis(R, D, d, t, SIGMA_O, E, mu=0.15)
    r_hi = jt.full_pullin_analysis(R, D, d, t, SIGMA_O, E, mu=0.35)
    assert r_hi['P1'] > r_lo['P1']


def test_tighter_bend_radius_gives_higher_plastic_strain():
    """A smaller R/d bends the riser harder -> more plastic strain -> k
    should increase as R/d decreases."""
    k_tight = jt.stage1_k(Rd=60.0, eps_o=0.002)
    k_loose = jt.stage1_k(Rd=200.0, eps_o=0.002)
    assert k_tight > k_loose
