# -*- coding: utf-8 -*-
"""
Reproduces the paper's design charts (Figs. 5-7) as curve FAMILIES
over a range of R/d, so the implementation is checked beyond the
single Appendix II point -- and gives the same design-chart tool the
paper's own "Design Application" section describes.

Saves PNGs to ../plots/ (relative to this script's own location, so
`python reproduce_figures.py` works standalone regardless of the
caller's working directory).

Fig. 5 (Stage I distances between contact points): L/d and l1/d vs
        R/d, one curve per D/d, at eps_o = 0.0025 (the paper's plotted
        material constant).
Fig. 6 (Critical J-tube bend angle): not reproduced -- the paper
        derives this from a SEPARATE stability criterion (a maximum
        pull-in load independent of bend angle beyond some critical
        value) that is not given as a closed-form equation in the
        body of the paper, only described qualitatively and shown as
        a chart. Flagged here rather than guessed.
Fig. 7 (Variation of pull-in load): (P1/Po)*1e3 vs R/d, one curve per
        D/d, at eps_o=0.0025, mu=0.3 -- Stage I friction+plasticity
        load only, matching the paper's own caption.
"""

import os

import numpy as np
import matplotlib.pyplot as plt

from walker1983_jtube import stage1_pullin

_OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'plots')
os.makedirs(_OUTDIR, exist_ok=True)

EPS_O = 0.0025          # plateau yield strain used throughout Figs. 5-7
SIGMA_O_E = EPS_O        # by definition, eps_o = sigma_o/E (paper's Nomenclature)


def _sweep(Dd_values, Rd_range, mu=None):
    """
    Sweeps R/d for each D/d in Dd_values, returning a dict
    D/d -> {'Rd':..., 'L_over_d':..., 'l1_over_d':..., 'P1_over_Po':...}
    skipping any R/d where the underlying formulas go non-physical
    (paper's own charts implicitly do the same -- they stop where the
    curves stop, they don't extrapolate through a singularity).
    """
    out = {}
    for Dd in Dd_values:
        Rds, Lds, l1ds, P1s = [], [], [], []
        for Rd in Rd_range:
            try:
                r = stage1_pullin(Rd, Dd, sigma_o=SIGMA_O_E, E=1.0,
                                   mu=(mu if mu is not None else 0.3))
                # sigma_o/E = SIGMA_O_E with E=1 is a trick to pass the
                # ratio directly without needing separate sigma_o, E
                # values -- stage1_pullin only ever uses their ratio.
            except ValueError:
                continue
            Rds.append(Rd)
            Lds.append(r['L_over_d'])
            l1ds.append(r['l1_over_d'])
            P1s.append(r['P1_over_Po'])
        out[Dd] = dict(Rd=np.array(Rds), L_over_d=np.array(Lds),
                        l1_over_d=np.array(l1ds), P1_over_Po=np.array(P1s))
    return out


def reproduce_fig5():
    Dd_values = [1.1, 1.3, 1.5, 1.7, 1.9, 2.1]
    Rd_range = np.linspace(100.0, 800.0, 141)
    sweep = _sweep(Dd_values, Rd_range)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 9))

    for Dd in Dd_values:
        s = sweep[Dd]
        if len(s['Rd']):
            ax1.plot(s['Rd'], s['L_over_d'], label=f'D/d={Dd}')
    ax1.set_xlabel('R/d')
    ax1.set_ylabel('L/d')
    ax1.set_title('Fig. 5 (upper) reproduction: Stage I -- L/d vs R/d\n'
                   '(eps_o = 0.0025)')
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.3)

    for Dd in Dd_values:
        s = sweep[Dd]
        if len(s['Rd']):
            ax2.plot(s['Rd'], s['l1_over_d'], label=f'D/d={Dd}')
    ax2.set_xlabel('R/d')
    ax2.set_ylabel('l1/d')
    ax2.set_title('Fig. 5 (lower) reproduction: Stage I -- l1/d vs R/d')
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    out = os.path.join(_OUTDIR, 'fig5_stage1_contact_distances.png')
    plt.savefig(out, dpi=150)
    plt.close(fig)
    print('Saved:', out)


def reproduce_fig7():
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 9))

    # upper panel: paper's Fig.7 upper -- D/d in {1.1..2.1}, mu=0.3
    Dd_values_upper = [1.1, 1.3, 1.5, 1.7, 1.9, 2.1]
    Rd_range_upper = np.linspace(100.0, 800.0, 141)
    sweep_upper = _sweep(Dd_values_upper, Rd_range_upper, mu=0.3)
    for Dd in Dd_values_upper:
        s = sweep_upper[Dd]
        if len(s['Rd']):
            ax1.plot(s['Rd'], s['P1_over_Po'] * 1e3, label=f'D/d={Dd}')
    ax1.set_xlabel('R/d')
    ax1.set_ylabel('(P1/Po) x 10^3')
    ax1.set_title('Fig. 7 (upper) reproduction: Stage I pull-in load\n'
                   '(eps_o=0.0025, mu=0.3)')
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.3)

    # lower panel: paper's Fig.7 lower -- D/d in {1.1, 1.5, 2.1}, R/d 0-200
    Dd_values_lower = [1.1, 1.5, 2.1]
    Rd_range_lower = np.linspace(20.0, 200.0, 91)
    sweep_lower = _sweep(Dd_values_lower, Rd_range_lower, mu=0.3)
    for Dd in Dd_values_lower:
        s = sweep_lower[Dd]
        if len(s['Rd']):
            ax2.plot(s['Rd'], s['P1_over_Po'] * 1e3, label=f'D/d={Dd}')
    ax2.set_xlabel('R/d')
    ax2.set_ylabel('(P1/Po) x 10^3')
    ax2.set_title('Fig. 7 (lower) reproduction: Stage I pull-in load, '
                   'lower R/d range')
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    out = os.path.join(_OUTDIR, 'fig7_stage1_pullin_load.png')
    plt.savefig(out, dpi=150)
    plt.close(fig)
    print('Saved:', out)


def reproduce_appendix2_marker():
    """
    Overlays the single Appendix II worked-example point onto a Fig.7-
    style curve, so the worked example's place on the design chart is
    visible directly rather than only as printed numbers.
    """
    Dd = 1.5
    Rd_range = np.linspace(20.0, 300.0, 141)
    sweep = _sweep([Dd], Rd_range, mu=0.3)[Dd]

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(sweep['Rd'], sweep['P1_over_Po'] * 1e3, color='#1a5c78',
            label=f'D/d={Dd}, mu=0.3 (this module)')

    r_ex = stage1_pullin(Rd=100.0, Dd=1.5, sigma_o=400e6, E=200e9, mu=0.3)
    ax.plot([100.0], [r_ex['P1_over_Po'] * 1e3], 'o', color='#97281f',
            markersize=9, zorder=5,
            label=f"Appendix II example (R/d=100): "
                  f"{r_ex['P1_over_Po']*1e3:.2f}")
    ax.axhline(2.2e-2 * 1e3, color='#97281f', linestyle=':', alpha=0.5,
               label='Paper-stated value: 22.0')

    ax.set_xlabel('R/d')
    ax.set_ylabel('(P1/Po) x 10^3')
    ax.set_title('Appendix II worked example on the Stage I design curve')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    out = os.path.join(_OUTDIR, 'appendix2_marker_on_design_curve.png')
    plt.savefig(out, dpi=150)
    plt.close(fig)
    print('Saved:', out)


if __name__ == '__main__':
    reproduce_fig5()
    reproduce_fig7()
    reproduce_appendix2_marker()
    print()
    print('Fig. 6 (critical J-tube bend angle) NOT reproduced: the paper')
    print('gives it only as a chart, derived from a stability criterion')
    print('not stated in closed form in the body of the paper. See this')
    print('script\'s module docstring.')
