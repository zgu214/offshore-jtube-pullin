# -*- coding: utf-8 -*-
"""
Walker, A.C., Davies, P. (1983). "A Design Basis for the J-Tube Method
of Riser Installation." ASME J. Energy Resources Technology, 105(3),
263-270. doi:10.1115/1.3264264

VERIFIED against the primary source (full paper PDF, read directly --
equation numbers below match the paper exactly). This implements the
paper's Stage I / II / III pull-in load analysis for a single-bend
J-tube (the "initial design" method the paper itself recommends -
see "Design Application").

-----------------------------------------------------------------------
Physical picture (see the paper's Fig. 2)
-----------------------------------------------------------------------
A pipeline riser is pulled through a curved J-tube of bend radius R,
inside diameter D, by a riser of outside diameter d. Three contact
points form a three-point bending condition; the riser between them is
treated as a built-in cantilever undergoing large elastic-plastic
deflection (bilinear material, Fig. 4). Three stages occur as the pull
proceeds:

  Stage I   - riser bends INTO the J-tube curvature. Governs initial
              design (paper: "considered sufficient... calculate the
              Stage I pull-in load only").
  Stage II  - pullhead is in the straight exit section; elastic
              relaxation reduces the load to a local minimum.
  Stage III - riser must straighten out its residual curvature; load
              rises again to about the same order as Stage I.

-----------------------------------------------------------------------
Equations implemented (numbers match the paper)
-----------------------------------------------------------------------
  (3)      M(omega) -- moment as plastic strain spreads from the
           outer fibres inward (quadrant model, Fig. 9)
  (4)      Mp = dm^2 * t * sigma_o ;  Mo = (pi/4) Mp
  (5),(6)  empirical fit to (3): plastic strain ratio eps/eps_o as a
           function of M/Mo (eq 5, valid eps/eps_o > 1) and the
           elastic branch (eq 6, eps/eps_o < 1)
  (7),(8)  alpha(k), beta(k) -- "plasticity correction factors"
           (k = Mbar/Mo), fitted to the elastic-plastic cantilever
           tip deflection/rotation (Fig. 12)
  (9),(10) cantilever tip rotation/deflection under elastic-plastic
           bending, in terms of alpha, beta
  (12)     L/d -- Stage I distance between the leading (B) and
           trailing (A) contact points; observed to stay CONSTANT
           through Stages I-III once established
  (13),(14) l1/d -- Stage I distance to the third (C) contact point
  (14')    the same, direct-entry case (paper's Fig. 13; no straight
           lead-in section) -- gives ~35% higher pull load in the
           worked example
  (15)     tau -- angle subtended at the bend centre by l1
  (2),(19) P1/Po -- Stage I pull-in load (plasticity + friction)
  (16),(17) Rp -- plastic radius of residual curvature after Stage I
  (20),(21) Stage II geometry: relative deflection delta, length l2,
           via constants C1, C2
  (22)-(24) Stage III: straightening deformation delta-hat, implicit
           alpha*Mbar/Mp vs l3/d relation, pull load P'/Po (chosen to
           MAXIMIZE P'/Po over l3 -- the riser "chooses" the length
           that most resists straightening)
  (25)     total Stage III load, P_III/Po = P_II/Po + P'/Po

Symbols follow the paper's Nomenclature. All lengths are carried as
dimensionless ratios (R/d, D/d, L/d, ...), matching how the paper's own
design charts (Figs. 5-7) are presented -- this also means the results
do not depend on absolute pipe size, only on the ratios.

-----------------------------------------------------------------------
Verification status
-----------------------------------------------------------------------
Cross-checked against the fully worked example in Appendix II (a
400 N/mm^2-yield riser, R/d=100, D/d=1.5, d=0.3 m, t=20 mm, mu=0.3):

    quantity        paper (graph-read)   this module (exact eq.)
    Po (kN)             7540                 7539.8  (<0.01%)
    Mbar/Mo (k)         --                   1.2351
    Mbar/Mp             0.97                 0.9701
    L/d                 32.7                 32.46   (0.7% low)
    l1/d                12.6                 12.68   (0.6% high)
    tau (deg)           7.1                  7.18    (1.1% high)
    P1/Po               2.2e-2               2.217e-2 (0.8% high)
    P1 (kN)             166                  167.2   (0.7% high)
    l2/d                34.5                 34.56   (0.2% high)
    P_II/Po             1.27e-2              1.292e-2 (1.7% high)
    Rp/d                166.7                166.67  (<0.1%)
    P_III/Po            2.19e-2              2.214e-2 (1.1% high)

Discrepancies of this size are expected and NOT bugs: the paper's own
worked example reads alpha*Mbar/Mp and beta*Mbar/Mp off Fig. 12 by eye
("from Fig. 12, alpha*Mbar/Mp = 1.11 and beta*Mbar/Mp = 1.07") rather
than evaluating eqs. (7)-(8) directly. This module evaluates the
equations directly, which is more precise than a hand graph-read, and
is the only way to compute results for parameter combinations the
paper didn't plot. See test_walker1983.py for the executable proof of
the table above (tolerances are set from this documented discrepancy,
not tightened to force a match).

ONE genuine inconsistency was found (not a graph-reading artifact):
the Fig. 7 caption, as transcribed from the scanned paper, defines
Po = pi*dm*t*sigma_o using the MEAN diameter dm = d-t, but the
Appendix II worked example's own arithmetic only reproduces the
stated 7540 kN using the OUTER diameter d (pi*dm*t*sigma_o gives
7037 kN, 6.7% low -- a real miss, not rounding). This module follows
d, because the worked example's number is independently verifiable
and a scanned figure-caption subscript is not. See
full_pullin_analysis()'s Po calculation.

Stage III's worked example is a coarser hand iteration (three trial
values of l3/d, picking the best by eye) and includes a borderline
elastic/plastic case (l3/d=17.5) that the paper treats as elastic by
inspection. This module instead solves alpha(k)*(k*pi/4) = target
uniformly via bisection for every trial, so results near that
boundary can differ from the paper's approximation by more than the
Stage I/II tolerance -- flagged explicitly in test_walker1983.py
rather than hidden behind a loose tolerance.

Not implemented:
  - multi-bend J-tubes (paper explicitly defers this to "a forthcoming
    paper" -- see the Concluding Remarks).
  - the back-tension/pull-cable-friction wrap
    P_T = (P1+T)*exp(mu_c*psi) + W (eq. after the Fig. 8 discussion) --
    a separate, simple multiplicative step the caller can apply on top
    of P1/P_T from here.
  - the direct-entry case (Fig. 13, no straight lead-in section ahead
    of the bend). full_pullin_analysis(direct_entry=True) raises
    NotImplementedError rather than guessing: the paper describes it
    in one paragraph ("l1 can be calculated using equation (14)") but
    Fig. 13's geometry visibly differs from the one eq (14) was
    derived for, and two straightforward reconstructions were tried
    and both missed the paper's own stated result (P1/Po=2.98e-2, a
    35% increase) by more than the ~1-2% every other number here
    achieves. See stage1_l1_over_d_direct_entry()'s docstring for what
    was tried.
"""

import math

import numpy as np


# =====================================================================
# eq (3)-(6): moment-curvature relationship
# =====================================================================

def moment_ratio_quadrant(omega):
    """
    eq (3): M/(dm^2 * t * sigma_o) = sin(omega) + sec(omega)*(pi/2 - omega)

    omega : half-angle (rad) subtended by the plastically-strained
            segments of the cross-section (Fig. 9). omega=0 -> first
            yield (M=Mo); omega=pi/2 -> fully plastic (M=Mp).

    Returns M / (dm^2 * t * sigma_o) -- multiply by dm^2*t*sigma_o for
    the moment in Nm (SI throughout).
    """
    omega = np.asarray(omega, dtype=float)
    out = np.empty_like(omega)
    at_limit = np.isclose(omega, math.pi / 2)
    out[at_limit] = 1.0                                    # sec(pi/2) sing.; limit = Mp
    ok = ~at_limit
    out[ok] = np.sin(omega[ok]) + (1.0 / np.cos(omega[ok])) * (math.pi / 2 - omega[ok])
    return out if out.shape else float(out)


def Mp_Mo(dm, t, sigma_o):
    """eq (4): fully-plastic and first-yield moments, Nm."""
    Mp = dm ** 2 * t * sigma_o
    Mo = (math.pi / 4) * Mp
    return Mp, Mo


def strain_ratio_from_moment(MMo, branch=None):
    """
    eq (5)/(6): (eps_bar/eps_o) as a function of (Mbar/Mo).

    eq (5) (plastic, eps/eps_o > 1):
        eps/eps_o = 1 + 2.25*(M/Mo - 1) + 4.42e5*(M/Mo - 1)**9
    eq (6) (elastic, eps/eps_o < 1):
        eps/eps_o = (4/pi)*(M/Mo)

    branch: force 'plastic' or 'elastic'; default auto-selects on
    MMo >= 1 (the two branches meet exactly at MMo=1, eps/eps_o=1).
    """
    if branch is None:
        branch = 'plastic' if MMo >= 1.0 else 'elastic'
    if branch == 'elastic':
        # eq (6) is stated in the paper as (4/pi)*(M/Mp), NOT (M/Mo).
        # Converting to this function's M/Mo argument: M/Mp = (M/Mo)*(Mo/Mp)
        # = MMo*(pi/4) [eq 4], so (4/pi)*(M/Mp) = (4/pi)*MMo*(pi/4) = MMo
        # exactly. I.e. expressed in M/Mo terms, eq (6) is just the trivial
        # linear-elastic identity eps/eps_o = M/Mo (strain proportional to
        # moment below first yield) -- which is also what
        # MMo_from_strain_ratio()'s elastic-branch shortcut already assumes.
        return MMo
    x = MMo - 1.0
    return 1.0 + 2.25 * x + 4.42e5 * x ** 9


def MMo_from_strain_ratio(eps_ratio, k_max=3.0, tol=1e-12):
    """
    Inverse of eq (5): given a target eps_bar/eps_o (> 1, plastic
    range), solve for k = Mbar/Mo by bisection. eq (5) is monotonic
    increasing in k for k >= 1, so bisection is unconditionally safe.
    """
    if eps_ratio <= 1.0:
        return eps_ratio  # elastic branch: eps/eps_o = M/Mo identically, see eq (6)
    lo, hi = 1.0, k_max
    f_hi = strain_ratio_from_moment(hi) - eps_ratio
    while f_hi < 0.0:
        hi *= 1.5
        f_hi = strain_ratio_from_moment(hi) - eps_ratio
        if hi > 1e6:
            raise ValueError("MMo_from_strain_ratio: no root found "
                              "(eps_ratio too large?)")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        f_mid = strain_ratio_from_moment(mid) - eps_ratio
        if abs(f_mid) < tol or (hi - lo) < tol:
            return mid
        if f_mid < 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# =====================================================================
# eq (7)-(8): plasticity correction factors alpha, beta
# =====================================================================

def alpha_factor(k):
    """
    eq (7): deflection plasticity correction factor.

        alpha = 2.25 - 1.875/k + 0.625/k^3
                + 1.31e5*(k-1)^10/k^3 + 1.21e5*(k-1)^11/k^3

    k = Mbar/Mo. For k <= 1 (no plastic straining at all) the cantilever
    is purely elastic and alpha = 1 EXACTLY by the definition in eq (1)
    (Delta = alpha*Delta_e) -- the polynomial above is only calibrated
    for k >= 1 (see Fig. 12's plotted domain) and is not evaluated
    below k=1.
    """
    k = np.asarray(k, dtype=float)
    out = np.where(
        k <= 1.0,
        1.0,
        2.25 - 1.875 / np.maximum(k, 1e-12) + 0.625 / np.maximum(k, 1e-12) ** 3
        + 1.31e5 * (k - 1.0) ** 10 / np.maximum(k, 1e-12) ** 3
        + 1.21e5 * (k - 1.0) ** 11 / np.maximum(k, 1e-12) ** 3,
    )
    return float(out) if out.shape == () else out


def beta_factor(k):
    """
    eq (8): rotation plasticity correction factor.

        beta = 2.25 - 2.5/k + 1.25/k^2 + 8.84e4*(k-1)^10/k^2

    Same k <= 1 -> beta = 1 convention as alpha_factor().
    """
    k = np.asarray(k, dtype=float)
    out = np.where(
        k <= 1.0,
        1.0,
        2.25 - 2.5 / np.maximum(k, 1e-12) + 1.25 / np.maximum(k, 1e-12) ** 2
        + 8.84e4 * (k - 1.0) ** 10 / np.maximum(k, 1e-12) ** 2,
    )
    return float(out) if out.shape == () else out


def k_from_alpha_Mbar_Mp(target, k_max=3.0, tol=1e-12):
    """
    Inverse problem used in Stage III: given a target value of the
    PRODUCT alpha(k)*(k*pi/4) [ = alpha*Mbar/Mp, the Fig. 12 y-axis
    quantity], solve for k = Mbar/Mo by bisection.

    alpha(k)*(k*pi/4) is monotonic increasing in k for k >= 1 (both
    factors increase), and equals pi/4 exactly at k=1 (alpha=1 there).
    For target <= pi/4 the section never left the elastic range, so
    alpha=1 identically and k = target*(4/pi) directly (no root-find
    needed, and none is well-posed since the product is constant at
    k<1 only in the trivial elastic sense that alpha=1 -- see
    alpha_factor's k<=1 convention).
    """
    if target <= math.pi / 4.0:
        return target * 4.0 / math.pi

    def g(k):
        return alpha_factor(k) * (k * math.pi / 4.0)

    lo, hi = 1.0, k_max
    g_hi = g(hi) - target
    while g_hi < 0.0:
        hi *= 1.5
        g_hi = g(hi) - target
        if hi > 1e6:
            raise ValueError("k_from_alpha_Mbar_Mp: no root found")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        gm = g(mid) - target
        if abs(gm) < tol or (hi - lo) < tol:
            return mid
        if gm < 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# =====================================================================
# Stage I geometry and load -- eq (12)-(15), (2)/(19)
# =====================================================================

def stage1_k(Rd, eps_o):
    """
    Stage I: the strain at contact point B is set purely by the
    J-tube's own curvature, eps_bar = d/(2R) (paper: "the strain at
    contact point B ... is the maximum strain"). Solve eq (5) for
    k = Mbar/Mo from that strain.

    Rd    : R/d (J-tube bend radius / riser OD)
    eps_o : plateau yield strain (sigma_o/E)
    """
    eps_ratio = (1.0 / (2.0 * Rd)) / eps_o
    return MMo_from_strain_ratio(eps_ratio)


def stage1_L_over_d(Rd, Dd, k, sigma_o_E):
    """
    eq (12): L/d, using alpha(k), beta(k) directly (NOT the Fig. 12
    product alpha*Mbar/Mp) together with Mbar/Mo = k:

        L/d = sqrt( (D/d - 1) /
                    [ (beta - 2*alpha/3)*k*(sigma_o/E)
                      - (beta*k)^2 * (sigma_o/E)^2 * (R/d + D/d)/2 ] )
    """
    a = alpha_factor(k)
    b = beta_factor(k)
    se = sigma_o_E
    denom = (b - 2.0 * a / 3.0) * k * se - (b * k) ** 2 * se ** 2 * (Rd + Dd) / 2.0
    if denom <= 0.0:
        raise ValueError(
            f"stage1_L_over_d: non-physical (denominator={denom:.3e} <= 0). "
            f"Check R/d, D/d, sigma_o/E are in a sensible range."
        )
    return math.sqrt((Dd - 1.0) / denom)


def stage1_l1_over_d(Rd, Dd, k, sigma_o_E):
    """
    eq (14): l1/d for the case with an adequate straight lead-in
    section ahead of the bend (the general/default case):

        l1/d = sqrt( [(R/d+D/d)^2 - (R/d+1)^2] /
                      [1 - (alpha*k)*(4*sigma_o/(3E))*(R/d+D/d)] )
    """
    a = alpha_factor(k)
    se = sigma_o_E
    num = (Rd + Dd) ** 2 - (Rd + 1.0) ** 2
    den = 1.0 - (a * k) * (4.0 * se / 3.0) * (Rd + Dd)
    if den <= 0.0:
        raise ValueError(f"stage1_l1_over_d: non-physical (denominator={den:.3e} <= 0)")
    return math.sqrt(num / den)


def stage1_l1_over_d_direct_entry(Rd, Dd, k, sigma_o_E):
    """
    NOT IMPLEMENTED -- see module docstring's "Not implemented" note.

    The paper covers the direct-entry case (Fig. 13: no straight
    lead-in section ahead of the bend) in one paragraph, stating only
    the headline RESULT for the Appendix II example (P1/Po = 2.98e-2,
    "an increase of 35 percent") and saying l1 "can be calculated
    using equation (14)" -- but Fig. 13's geometry visibly differs
    from Fig. 2(a)'s three-point layout that eq (14) was derived for
    (Fig. 13 only shows two dimensions, l1 and L, not three contact
    points), so eq (14) cannot be the WHOLE story.

    Two straightforward readings were tried and BOTH rejected because
    neither reproduces 2.98e-2 to the ~1-2% tolerance every other
    number in this module hits:

      (a) reuse eq (19) with L replaced by l1 (no separate lead-in
          length to establish L) -> gives 3.13e-2, 5% high
      (b) keep L/d=32.7 as in the standard case, solve backwards for
          whatever l1/d would be needed -> requires l1/d=8.3, which is
          not derivable from eq (14) or anything else stated

    Guessing a third formula to force a match to one worked number
    would be fitting, not verification. Raising here rather than
    silently returning stage1_l1_over_d()'s unmodified value (which
    would give an IDENTICAL result to the lead-in case -- silently
    wrong, since the paper is explicit the two cases differ by ~35%).
    """
    raise NotImplementedError(
        "Direct-entry J-tube geometry (Fig. 13) is not implemented: the "
        "source paper's description ('l1 can be calculated using equation "
        "(14)') under-specifies a geometry that visibly differs from the "
        "one eq (14) was derived for, and no reconstruction tried "
        "reproduces the paper's own stated result (P1/Po=2.98e-2) to "
        "within the tolerance every other quantity in this module "
        "achieves. See this function's docstring for what was tried. "
        "Only the standard lead-in-section case (Fig. 2a, the default) "
        "is implemented."
    )


def stage1_tau_deg(Rd, Dd, l1_over_d):
    """eq (15): angle (deg) subtended by l1 at the bend centre O."""
    ratio = l1_over_d / (Rd + Dd)
    return math.degrees(math.asin(min(1.0, ratio)))


def pullin_load_ratio(Mbar_Mp, Rd, sigma_o_E, mu, L_over_d, K_over_d):
    """
    eq (2)/(19): general nondimensional pull-in load,

        P/Po = (Mbar/(pi*Mp)) * [ d/R - 2*sigma_o/E + 2*mu*(d/L + K) ]

    The (Mbar/pi/Mp) prefactor multiplies the WHOLE bracket, including
    the friction term -- this is easy to mis-parse from the paper's
    prose; the Appendix II worked example is only reproduced correctly
    with the prefactor applied to both terms (see module docstring).

    Rd       : R/d
    L_over_d : the (stage-independent, fixed-after-Stage-I) L/d
    K_over_d : the stage-specific "K" term, already in the same
               d/(...) form as the paper's K = d/l1 / d/l2 / etc.
               (i.e. pass K_over_d = 1/(l1/d) for Stage I, etc.)
    """
    bracket = (1.0 / Rd) - 2.0 * sigma_o_E + 2.0 * mu * (1.0 / L_over_d + K_over_d)
    return (Mbar_Mp / math.pi) * bracket


def stage1_pullin(Rd, Dd, sigma_o, E, mu, direct_entry=False):
    """
    Full Stage I solve: strain -> k -> L/d, l1/d, tau -> P1/Po.

    Returns a dict with every intermediate quantity (so a caller can
    cross-check against the paper's worked example term by term, as
    test_walker1983.py does).
    """
    sigma_o_E = sigma_o / E
    eps_o = sigma_o_E  # plateau yield strain = plateau yield stress / E (paper's def.)

    k = stage1_k(Rd, eps_o)
    Mbar_Mp = k * math.pi / 4.0
    a = alpha_factor(k)
    b = beta_factor(k)

    Ld = stage1_L_over_d(Rd, Dd, k, sigma_o_E)
    if direct_entry:
        l1d = stage1_l1_over_d_direct_entry(Rd, Dd, k, sigma_o_E)
    else:
        l1d = stage1_l1_over_d(Rd, Dd, k, sigma_o_E)
    tau = stage1_tau_deg(Rd, Dd, l1d)

    P1_Po = pullin_load_ratio(Mbar_Mp, Rd, sigma_o_E, mu, Ld, 1.0 / l1d)

    return dict(
        k=k, Mbar_Mo=k, Mbar_Mp=Mbar_Mp, alpha=a, beta=b,
        L_over_d=Ld, l1_over_d=l1d, tau_deg=tau,
        P1_over_Po=P1_Po, sigma_o_E=sigma_o_E,
    )


# =====================================================================
# Stage II -- eq (20), (21)
# =====================================================================

def stage2_C1_C2(Rd, sigma_o_E, no_plastic_straining=False):
    """
    C1, C2 as defined beneath eq (21).

        general:              C1 = (d/R - sigma_o/E)^-1
                               C2 = (d/R - 2*sigma_o/(3E))^-1
        no plastic straining: C1 = R/(2d) ; C2 = 2R/(3d)
    """
    if no_plastic_straining:
        return Rd / 2.0, 2.0 * Rd / 3.0
    C1 = 1.0 / (1.0 / Rd - sigma_o_E)
    C2 = 1.0 / (1.0 / Rd - 2.0 * sigma_o_E / 3.0)
    return C1, C2


def stage2_l2_over_d(Rd, Dd, C1, C2):
    """
    eq (21): l2/d between the leading contact points as the pullhead
    exits into the straight section.

        l2/d = sqrt( 2*C1*(D/d-1) / [(2 - C1/C2) - (R/d+D/d)/C1] )
    """
    denom = (2.0 - C1 / C2) - (Rd + Dd) / C1
    if denom <= 0.0:
        raise ValueError(f"stage2_l2_over_d: non-physical (denominator={denom:.3e} <= 0)")
    return math.sqrt(2.0 * C1 * (Dd - 1.0) / denom)


def stage2_pullin(Rd, Dd, sigma_o, E, mu, k, L_over_d):
    """
    Full Stage II solve. k (=Mbar/Mo) and L/d are carried over
    unchanged from Stage I -- the paper's tests show L stays constant
    through the pull, and Mbar depends only on the (unchanged)
    J-tube bend radius, not on which stage the pullhead is in.
    """
    sigma_o_E = sigma_o / E
    Mbar_Mp = k * math.pi / 4.0

    C1, C2 = stage2_C1_C2(Rd, sigma_o_E)
    l2d = stage2_l2_over_d(Rd, Dd, C1, C2)

    P2_Po = pullin_load_ratio(Mbar_Mp, Rd, sigma_o_E, mu, L_over_d, 1.0 / l2d)

    return dict(C1=C1, C2=C2, l2_over_d=l2d, P2_over_Po=P2_Po)


# =====================================================================
# Stage III -- eq (16), (17), (22)-(25)
# =====================================================================

def plastic_radius_of_curvature(Rd, sigma_o_E):
    """
    eq (17): d/Rp = d/R - 2*sigma_o/E  ->  Rp/d.

    Rp is the residual (permanently bent) radius of curvature left in
    the riser after Stage I plastic straining and elastic spring-back.
    """
    inv_Rpd = (1.0 / Rd) - 2.0 * sigma_o_E
    if inv_Rpd <= 0.0:
        raise ValueError("plastic_radius_of_curvature: riser did not "
                          "strain plastically enough for a residual "
                          "curvature to exist (d/Rp <= 0)")
    return 1.0 / inv_Rpd


def stage3_alpha_Mbar_Mp(l3_over_d, Rp_over_d, Dd, sigma_o_E):
    """
    eq (22)+(23) combined into a single closed form for
    alpha*Mbar/Mp as a function of a TRIAL l3/d.

    eq (22): delta_hat = l3^2/(2*Rp) - (D - d)
    eq (23): delta_hat/d = (alpha*Mbar/Mp)*(8*sigma_o/(3*pi*E))*(l3/d)^2

    Combining and solving for (alpha*Mbar/Mp):

        alpha*Mbar/Mp = [ (l3/d)^2/(2*Rp/d) - (D/d - 1) ]
                        / [ (8*sigma_o/(3*pi*E)) * (l3/d)^2 ]
    """
    num = (l3_over_d ** 2) / (2.0 * Rp_over_d) - (Dd - 1.0)
    den = (8.0 * sigma_o_E / (3.0 * math.pi)) * l3_over_d ** 2
    return num / den


def stage3_P_prime_over_Po(Mbar_Mp, mu, l3_over_d):
    """eq (24): P' (straightening load) / Po = (2*mu/pi)*(Mbar/Mp)*(d/l3)."""
    return (2.0 * mu / math.pi) * Mbar_Mp * (1.0 / l3_over_d)


def stage3_pullin(Rd, Dd, sigma_o, E, mu, k, P2_over_Po,
                   l3_over_d_trials=None):
    """
    Full Stage III solve: for each trial l3/d, get alpha*Mbar/Mp from
    the closed form (eq 22-23), invert it to Mbar/Mp (eq 7's alpha(k)
    relation), then compute P'/Po (eq 24). The riser "selects" l3 to
    MAXIMIZE P'/Po (paper: "the length l3 ... must be chosen to
    maximize the value of P'/Po") -- so we scan and report the peak.

    Default trial range follows the worked example's exploratory
    grid (roughly d/R-scaled), refined to a finer scan so the peak is
    actually located rather than eyeballed from 3 points as the paper
    does by hand.

    Returns the peak-P' result plus the full scan (for plotting/audit).
    P_III/Po = P2_over_Po + max(P'/Po), per eq (25).
    """
    sigma_o_E = sigma_o / E
    Rpd = plastic_radius_of_curvature(Rd, sigma_o_E)

    if l3_over_d_trials is None:
        l3_over_d_trials = np.linspace(5.0, 60.0, 221)

    rows = []
    for l3d in l3_over_d_trials:
        target = stage3_alpha_Mbar_Mp(l3d, Rpd, Dd, sigma_o_E)
        if target <= 0.0:
            continue  # geometrically inadmissible (delta_hat <= 0): no straightening work yet
        k3 = k_from_alpha_Mbar_Mp(target)
        Mbar_Mp3 = k3 * math.pi / 4.0
        Pprime_Po = stage3_P_prime_over_Po(Mbar_Mp3, mu, l3d)
        rows.append((l3d, target, k3, Mbar_Mp3, Pprime_Po))

    if not rows:
        raise ValueError("stage3_pullin: no admissible l3/d in the trial range")

    arr = np.array(rows)
    i_peak = int(np.argmax(arr[:, 4]))
    l3d_peak, target_peak, k3_peak, Mbar_Mp3_peak, Pprime_Po_peak = arr[i_peak]

    return dict(
        Rp_over_d=Rpd,
        l3_over_d=l3d_peak, alpha_Mbar_Mp=target_peak,
        k=k3_peak, Mbar_Mp=Mbar_Mp3_peak,
        P_prime_over_Po=Pprime_Po_peak,
        P3_over_Po=P2_over_Po + Pprime_Po_peak,
        scan=arr,  # columns: l3/d, alpha*Mbar/Mp, k, Mbar/Mp, P'/Po
    )


# =====================================================================
# Convenience: full three-stage solve from raw engineering inputs
# =====================================================================

def full_pullin_analysis(R, D, d, t, sigma_o, E, mu, direct_entry=False):
    """
    Top-level entry point taking ENGINEERING units (m, Pa) rather than
    the dimensionless ratios the low-level functions use.

    R        : J-tube bend radius (m)
    D        : J-tube inside diameter (m)
    d        : riser outside diameter (m)
    t        : riser wall thickness (m)
    sigma_o  : riser plateau yield stress (Pa)
    E        : riser Young's modulus (Pa)
    mu       : riser/J-tube coefficient of friction

    Returns a dict with 'stage1', 'stage2', 'stage3' sub-dicts (ratios)
    plus 'Po' (N) and each stage's dimensional pull load (N).
    """
    Rd, Dd = R / d, D / d
    # Po = pi*d*t*sigma_o, using the OUTER diameter d.
    #
    # The Fig. 7 caption (as transcribed from the scanned paper) reads
    # "Po = pi*dm*t*sigma_o" with the MEAN diameter dm = d-t. But the
    # Appendix II worked example's own arithmetic only reproduces
    # exactly (7540 kN, matched to the kN) using the OUTER diameter d:
    #     pi*d*t*sigma_o  = pi*0.30*0.02*400e6 = 7539.8 kN  (paper: 7540)
    #     pi*dm*t*sigma_o = pi*0.28*0.02*400e6 = 7037.2 kN  (6.7% low)
    # The worked example is verifiable ground truth; a scanned
    # subscript is not, so this module follows d. If your own project
    # data was fitted/calibrated against the dm convention instead,
    # override by passing d - t here.
    Po = math.pi * d * t * sigma_o

    s1 = stage1_pullin(Rd, Dd, sigma_o, E, mu, direct_entry=direct_entry)
    s2 = stage2_pullin(Rd, Dd, sigma_o, E, mu, s1['k'], s1['L_over_d'])
    s3 = stage3_pullin(Rd, Dd, sigma_o, E, mu, s1['k'], s2['P2_over_Po'])

    return dict(
        Rd=Rd, Dd=Dd, Po=Po,
        stage1=s1, stage2=s2, stage3=s3,
        P1=s1['P1_over_Po'] * Po,
        P2=s2['P2_over_Po'] * Po,
        P3=s3['P3_over_Po'] * Po,
    )


# =====================================================================
# Standalone demo: reproduce the Appendix II worked example
# =====================================================================

if __name__ == '__main__':
    # Appendix II inputs, exactly as stated in the paper.
    sigma_o = 400e6      # Pa (400 N/mm^2)
    E = 200e9            # Pa
    R = 30.0             # m  (R/d = 100)
    D = 0.45             # m  (D/d = 1.5)
    d = 0.30             # m
    t = 0.020            # m
    mu = 0.30

    result = full_pullin_analysis(R, D, d, t, sigma_o, E, mu)

    print('=' * 72)
    print('Walker & Davies (1983) -- Appendix II worked example, reproduced')
    print('=' * 72)
    print(f'  Po           = {result["Po"]/1e3:8.1f} kN   (paper: 7540 kN)')
    print()
    s1 = result['stage1']
    print('  --- Stage I ---')
    print(f'  Mbar/Mo (k)  = {s1["k"]:8.4f}')
    print(f'  Mbar/Mp      = {s1["Mbar_Mp"]:8.4f}   (paper Fig.10: 0.97)')
    print(f'  L/d          = {s1["L_over_d"]:8.2f}   (paper: 32.7)')
    print(f'  l1/d         = {s1["l1_over_d"]:8.2f}   (paper: 12.6)')
    print(f'  tau (deg)    = {s1["tau_deg"]:8.2f}   (paper: 7.1)')
    print(f'  P1/Po        = {s1["P1_over_Po"]:.4e} (paper: 2.2e-2)')
    print(f'  P1           = {result["P1"]/1e3:8.1f} kN   (paper: 166 kN)')
    print()
    s2 = result['stage2']
    print('  --- Stage II ---')
    print(f'  C1, C2       = {s2["C1"]:8.2f}, {s2["C2"]:8.2f}   (paper: 125, 115.4)')
    print(f'  l2/d         = {s2["l2_over_d"]:8.2f}   (paper: 34.5)')
    print(f'  P2/Po        = {s2["P2_over_Po"]:.4e} (paper: 1.27e-2)')
    print()
    s3 = result['stage3']
    print('  --- Stage III ---')
    print(f'  Rp/d         = {s3["Rp_over_d"]:8.2f}   (paper: 166.7)')
    print(f'  l3/d (peak)  = {s3["l3_over_d"]:8.2f}   (paper explored 17.5/20/22.5)')
    print(f"  P'/Po (peak) = {s3['P_prime_over_Po']:.4e} (paper best trial: 0.92e-2)")
    print(f'  P_III/Po     = {s3["P3_over_Po"]:.4e} (paper: 2.19e-2)')
    print(f'  P_III        = {result["P3"]/1e3:8.1f} kN')
    print('=' * 72)
    print('Discrepancies of order 1-2% vs. the paper are expected: the ')
    print('paper reads alpha*Mbar/Mp, beta*Mbar/Mp off a hand-drawn chart')
    print('(Fig. 12); this evaluates eqs. (7)-(8) directly. See the module')
    print('docstring and test_walker1983.py.')
