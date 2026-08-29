# -*- coding: utf-8 -*-
"""
Tests for regime_classifier.py.

Run with:  pytest test_regime_classifier.py -v
"""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pytest

import regime_classifier as rc


# =====================================================================
# hollow_tube_second_moment -- standard mechanics-of-materials result
# =====================================================================

def test_hollow_tube_second_moment_matches_closed_form():
    OD, t = 0.30, 0.020
    ID = OD - 2 * t
    expected = (math.pi / 64.0) * (OD ** 4 - ID ** 4)
    assert rc.hollow_tube_second_moment(OD, t) == pytest.approx(expected, rel=1e-12)


def test_hollow_tube_second_moment_solid_rod_limit():
    """t = OD/2 (just short of it) -> ID -> 0 -> a solid rod, I -> pi/64*OD^4."""
    OD = 0.30
    t = OD / 2.0 - 1e-9
    I = rc.hollow_tube_second_moment(OD, t)
    assert I == pytest.approx((math.pi / 64.0) * OD ** 4, rel=1e-6)


def test_hollow_tube_second_moment_rejects_impossible_geometry():
    with pytest.raises(ValueError):
        rc.hollow_tube_second_moment(OD=0.30, t=0.20)   # t > OD/2
    with pytest.raises(ValueError):
        rc.hollow_tube_second_moment(OD=0.30, t=0.0)    # zero wall
    with pytest.raises(ValueError):
        rc.hollow_tube_second_moment(OD=0.30, t=-0.01)  # negative wall


def test_hollow_tube_second_moment_increases_with_wall_thickness():
    OD = 0.30
    I_thin = rc.hollow_tube_second_moment(OD, t=0.005)
    I_thick = rc.hollow_tube_second_moment(OD, t=0.020)
    assert I_thick > I_thin


# =====================================================================
# wrap_threshold_tension
# =====================================================================

def test_wrap_threshold_matches_manual_EI_computation():
    """Confirms the convenience wrapper computes EI = E*I exactly as a
    caller doing it by hand would, not some other combination."""
    OD, t, E, r, D = 0.30, 0.020, 200e9, 30.0, 0.45
    EI_manual = E * rc.hollow_tube_second_moment(OD, t)
    a = None
    import ronalds1992_jtube as rn
    a = rn.straight_approach_length(r, D, OD)
    expected = rn.can_umbilical_wrap_bend(EI_manual, a, r, D, OD, f=4.0)
    assert rc.wrap_threshold_tension(OD, t, E, r, D) == pytest.approx(expected, rel=1e-12)


def test_wrap_threshold_increases_with_stiffer_wall():
    """A stiffer (thicker-walled) pipe needs MORE tension to be forced
    into wrapping the bend -- basic physical direction."""
    r, D, E = 30.0, 0.45, 200e9
    T_thin = rc.wrap_threshold_tension(OD=0.30, t=0.005, E=E, r=r, D=D)
    T_thick = rc.wrap_threshold_tension(OD=0.30, t=0.020, E=E, r=r, D=D)
    assert T_thick > T_thin


def test_wrap_threshold_decreases_with_larger_bend_radius():
    """A larger (gentler) bend radius makes wrapping easier -- lower
    tension needed -- since curvature develops more gradually."""
    OD, t, E, D = 0.30, 0.020, 200e9, 0.45
    T_tight = rc.wrap_threshold_tension(OD, t, E, r=12.0, D=D)
    T_gentle = rc.wrap_threshold_tension(OD, t, E, r=60.0, D=D)
    assert T_gentle < T_tight


# =====================================================================
# classify_pullin_regime -- the headline function
# =====================================================================

def test_classify_without_pullin_tension_gives_no_regime_label():
    """A threshold alone doesn't tell you which side of it you're on --
    must not silently guess a regime when no tension is supplied."""
    r = rc.classify_pullin_regime(OD=0.30, t=0.020, E=200e9, r=30.0, D=0.45)
    assert r['regime'] == 'unknown (no pull tension given)'
    assert r['margin'] is None
    assert r['threshold_T'] > 0


def test_classify_rejects_nonpositive_pullin_tension():
    with pytest.raises(ValueError):
        rc.classify_pullin_regime(OD=0.30, t=0.020, E=200e9, r=30.0, D=0.45,
                                  pullin_tension=0.0)
    with pytest.raises(ValueError):
        rc.classify_pullin_regime(OD=0.30, t=0.020, E=200e9, r=30.0, D=0.45,
                                  pullin_tension=-1000.0)


def test_walker_appendix_ii_riser_classifies_as_rigid():
    """
    THE key check: Walker & Davies' own worked-example riser, at its
    own computed Stage I pull load, should classify unambiguously as
    'rigid (Walker-type)' -- confirming the paper's own example is
    self-consistent with the discrete-contact assumption it relies on,
    not accidentally sitting in Ronalds' regime instead.
    """
    r = rc.classify_pullin_regime(
        OD=0.30, t=0.020, E=200e9, r=30.0, D=0.45,
        pullin_tension=167.2e3,   # Walker's own computed P1 (see walker1983_jtube.py)
    )
    assert r['margin'] == pytest.approx(11.5, rel=0.05)
    assert 'rigid' in r['regime']


def test_a_much_more_flexible_object_classifies_as_flexible():
    """
    A placeholder umbilical-scale bending stiffness (orders of
    magnitude below Walker's steel riser), pulled at a typical
    umbilical-scale tension, should classify as 'flexible
    (Ronalds-type)' -- confirming the classifier actually discriminates
    both ways, not just always saying 'rigid'.
    """
    r = rc.classify_pullin_regime(
        OD=0.10, t=0.010, E=1e9,   # a soft, thin-wall placeholder -- NOT steel
        r=12.0, D=0.20,
        pullin_tension=50e3,        # 50 kN, a plausible umbilical pull tension
    )
    assert 'flexible' in r['regime']


def test_marginal_case_between_the_two_regimes():
    """Construct a case right at the threshold (margin ~= 1) and
    confirm it lands in the explicitly-uncertain 'marginal' bucket
    rather than confidently asserting either regime."""
    r_geom, D = 30.0, 0.45
    # find an OD/t/E combination whose threshold is close to a chosen
    # tension by using the same placeholder object but tuning tension
    # to sit right at its own threshold
    probe = rc.classify_pullin_regime(OD=0.30, t=0.020, E=200e9,
                                      r=r_geom, D=D, pullin_tension=1.0)
    threshold = probe['threshold_T']

    at_threshold = rc.classify_pullin_regime(
        OD=0.30, t=0.020, E=200e9, r=r_geom, D=D, pullin_tension=threshold)
    assert at_threshold['margin'] == pytest.approx(1.0, rel=1e-9)
    assert 'marginal' in at_threshold['regime']


def test_classify_result_dict_has_expected_keys():
    r = rc.classify_pullin_regime(OD=0.30, t=0.020, E=200e9, r=30.0, D=0.45,
                                  pullin_tension=100e3)
    assert set(r.keys()) == {'EI', 'threshold_T', 'pullin_tension', 'margin', 'regime'}
