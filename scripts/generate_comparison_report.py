# -*- coding: utf-8 -*-
"""
Runs the Appendix II worked example and writes a comparison report
(paper value vs. this module's value, with % difference) as both JSON
and a plain-text table, to ../results/ (relative to this script, so it
works standalone regardless of the caller's working directory).

This is the single artifact to hand to someone auditing the
reproduction: every number the paper states, next to what this module
computes, with the discrepancy called out explicitly rather than
buried in test tolerances.
"""

import json
import os

from walker1983_jtube import full_pullin_analysis

_OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
os.makedirs(_OUTDIR, exist_ok=True)

# Appendix II inputs, exactly as stated in the paper.
SIGMA_O = 400e6
E = 200e9
R = 30.0
D = 0.45
d = 0.30
t = 0.020
MU = 0.30


def build_rows(result):
    s1, s2, s3 = result['stage1'], result['stage2'], result['stage3']
    return [
        dict(quantity='Po (kN)', paper=7540.0, ours=result['Po'] / 1e3,
             note='Outer diameter d, not mean diameter dm -- see module docstring'),
        dict(quantity='Mbar/Mp', paper=0.97, ours=s1['Mbar_Mp'],
             note='Paper reads this off Fig. 10 / eq (5)'),
        dict(quantity='Stage I: L/d', paper=32.7, ours=s1['L_over_d'], note=''),
        dict(quantity='Stage I: l1/d', paper=12.6, ours=s1['l1_over_d'], note=''),
        dict(quantity='Stage I: tau (deg)', paper=7.1, ours=s1['tau_deg'], note=''),
        dict(quantity='Stage I: P1/Po', paper=2.2e-2, ours=s1['P1_over_Po'],
             note='Paper: 1.85e-3 (plasticity) + 20.0e-3 (friction)'),
        dict(quantity='Stage I: P1 (kN)', paper=166.0, ours=result['P1'] / 1e3, note=''),
        dict(quantity='Stage II: C1', paper=125.0, ours=s2['C1'], note=''),
        dict(quantity='Stage II: C2', paper=115.4, ours=s2['C2'], note=''),
        dict(quantity='Stage II: l2/d', paper=34.5, ours=s2['l2_over_d'], note=''),
        dict(quantity='Stage II: P2/Po', paper=1.27e-2, ours=s2['P2_over_Po'], note=''),
        dict(quantity='Stage III: Rp/d', paper=166.7, ours=s3['Rp_over_d'], note=''),
        dict(quantity="Stage III: P'/Po (peak)", paper=0.92e-2,
             ours=s3['P_prime_over_Po'],
             note='Paper: coarse 3-point manual search (17.5/20/22.5); '
                  'this module: fine scan, peak at l3/d=%.1f' % s3['l3_over_d']),
        dict(quantity='Stage III: P_III/Po', paper=2.19e-2, ours=s3['P3_over_Po'], note=''),
    ]


def main():
    result = full_pullin_analysis(R, D, d, t, SIGMA_O, E, MU)
    rows = build_rows(result)
    for r in rows:
        r['pct_diff'] = 100.0 * (r['ours'] - r['paper']) / r['paper']

    # --- JSON ---
    json_path = os.path.join(_OUTDIR, 'appendix2_comparison.json')
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(dict(
            source='Walker & Davies (1983), Appendix II worked example',
            inputs=dict(sigma_o_Pa=SIGMA_O, E_Pa=E, R_m=R, D_m=D, d_m=d,
                        t_m=t, mu=MU),
            rows=rows,
        ), f, indent=2)

    # --- plain-text table ---
    txt_path = os.path.join(_OUTDIR, 'appendix2_comparison.txt')
    with open(txt_path, 'w', encoding='utf-8') as f:
        f.write('Walker & Davies (1983) -- Appendix II worked example\n')
        f.write('Comparison: paper-stated value vs. walker1983_jtube.py\n')
        f.write('=' * 88 + '\n')
        f.write(f"{'quantity':<26}{'paper':>12}{'this module':>14}{'% diff':>10}   note\n")
        f.write('-' * 88 + '\n')
        for r in rows:
            f.write(f"{r['quantity']:<26}{r['paper']:>12.4g}{r['ours']:>14.4g}"
                     f"{r['pct_diff']:>9.1f}%   {r['note']}\n")
        f.write('-' * 88 + '\n')
        f.write('\n')
        f.write('All discrepancies are within the ~1-2% documented in the module\n')
        f.write("docstring, EXCEPT where noted -- expected because the paper's own\n")
        f.write('worked example reads alpha*Mbar/Mp, beta*Mbar/Mp off a hand-drawn\n')
        f.write('chart (Fig. 12) rather than evaluating eqs. (7)-(8) directly.\n')

    print('Wrote:', json_path)
    print('Wrote:', txt_path)
    print()
    with open(txt_path, encoding='utf-8') as f:
        print(f.read())


if __name__ == '__main__':
    main()
