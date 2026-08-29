# J-tube pull-in analysis

Python reproductions of two J-tube pull-in design methods:

- **Walker, A.C., Davies, P. (1983)**, "A Design Basis for the J-Tube
  Method of Riser Installation," *ASME J. Energy Resources Technology*,
  105(3), 263-270 -- elastic-plastic RIGID PIPE bending, single bend.
- **Ronalds, B.F. (1992)**, "J-Tube Design for Flexible Umbilicals,"
  OTC 6875 -- cable/friction theory for a FLEXIBLE umbilical, single
  and multi-bend. Explicitly builds on Walker & Davies (cites it as
  Ref. [1]) and extends to the multi-bend case Walker & Davies defers.

## Are the two theories the same? No -- and that's checkable, not just asserted

The two papers make **opposite assumptions about how the object
contacts the J-tube wall**, and that single difference drives
everything else:

| | Walker & Davies | Ronalds |
|---|---|---|
| Contact | Discrete -- 3 points, straight between them | Continuous -- hugs the wall throughout the bend |
| Load mechanism | Elastic-plastic bending (dominant) + friction (secondary) | Pure friction (capstan equation); bending negligible by assumption |
| `sigma_o`, `E` | Central -- set how much the pipe plastically deforms | Peripheral -- only used in a secondary "can it even wrap?" check |
| Valid for | Stiff risers (`d_m/t < 30`, explicitly does NOT conform to the bend) | Flexible umbilicals (bending stiffness assumed negligible) |

**You cannot use one paper's numbers to validate the other's method in
general** -- they compute different physical quantities under
different governing assumptions, and plugging a stiff riser's
geometry into the flexible-cable formulas just answers a different
(and, for that riser, wrong) question.

What genuinely IS shared: `R/d`, `D/d`, bend angle, and `mu` carry the
same physical meaning in both. And Ronalds' pure capstan equation is
the correct limit of Walker's friction term as bending stiffness -> 0
-- a real theoretical link, but a limit, not an identity over the
range either paper is meant for. Their domains of applicability are
essentially disjoint, not overlapping.

`regime_classifier.py` makes this checkable rather than asserted: it
computes the tension above which a given object (OD, wall thickness,
E) would actually conform continuously to a given J-tube bend
(Ronalds' own wrap-vs-bridge criterion, eq. 15), and compares it
against an estimated pull tension. Run against Walker's own Appendix
II riser at its own computed Stage I load, the wrap threshold is
**~11.5x higher** than the tension that riser is actually pulled with
-- confirming the paper's own worked example is self-consistent with
the discrete-contact assumption it depends on, nowhere near Ronalds'
regime.

### Combined GUI

`combined_gui.py` puts all three tools -- Walker, Ronalds, and the
Regime Classifier -- in one window as three tabs, so a single geometry
can flow between them without retyping numbers. It does not change any
calculation logic from the standalone GUIs (which remain fully usable
on their own); it only re-hosts the same input/compute/display pattern
inside shared tabs. The one new piece of behaviour: the Regime
Classifier tab has a **"Use Walker tab's P1"** button that pulls the
currently-computed Stage I pull load straight from the Walker tab into
the classifier's pull-tension field -- the same manual step that
produced the ~11.5x margin finding above, now one click instead of
copying a number by hand.

## Folder layout

```
J_tube_2026/
├── scripts/
│   ├── walker1983_jtube.py            core library, Walker & Davies (1983)
│   ├── test_walker1983.py             30 tests
│   ├── walker1983_gui.py              desktop GUI (tkinter + matplotlib)
│   ├── ronalds1992_jtube.py           core library, Ronalds (1992)
│   ├── test_ronalds1992.py            28 tests
│   ├── ronalds1992_gui.py             desktop GUI, two tabs (multi-bend / single-bend support)
│   ├── reproduce_figures.py           regenerates Walker Figs. 5 and 7
│   ├── generate_comparison_report.py  writes the Walker paper-vs-module table
│   ├── regime_classifier.py           cross-paper: which paper's assumptions actually apply
│   ├── test_regime_classifier.py      13 tests
│   └── combined_gui.py                all three tools in one window (see below)
├── results/
│   ├── appendix2_comparison.json      machine-readable comparison (Walker)
│   └── appendix2_comparison.txt       the same, as a readable table
├── plots/
│   ├── fig5_stage1_contact_distances.png
│   ├── fig7_stage1_pullin_load.png
│   └── appendix2_marker_on_design_curve.png
└── *.pdf                              source papers (see "Other papers" below)
```

## Quick start

```bash
cd scripts
python -m pytest test_walker1983.py test_ronalds1992.py test_regime_classifier.py -v   # 71 tests, all pass
python walker1983_jtube.py                    # prints the Appendix II worked example
python ronalds1992_jtube.py                    # prints the multi-bend/guide-support demo
python regime_classifier.py                     # which paper's assumptions apply, for a given object
python walker1983_gui.py                        # desktop GUI, pre-filled with Appendix II
python ronalds1992_gui.py                        # desktop GUI, multi-bend + single-bend support tabs
python combined_gui.py                          # all three tools in one window (recommended)
python generate_comparison_report.py             # writes results/appendix2_comparison.*
python reproduce_figures.py                      # writes plots/*.png
```

Every script resolves its own output path relative to its own file
location, so all of these run correctly regardless of your current
working directory.

## Walker & Davies (1983) -- what's implemented

The paper's full single-bend, three-stage pull-in analysis:

| Stage | What it models | Status |
|---|---|---|
| I | Riser bending into the J-tube curvature (governs initial design) | Verified, ~0.7% of paper |
| II | Elastic relaxation as the pullhead enters the straight exit | Verified, ~1.7% of paper |
| III | Riser straightening its residual curvature | Verified, ~1.1% of paper (peak found by fine scan, not the paper's 3-point hand search) |

Plus: the moment-curvature relationship (eq. 3-6), the elastic-plastic
"plasticity correction factors" alpha/beta (eq. 7-8, Fig. 12), and
design-chart reproduction of Figs. 5 and 7.

**Not implemented**, both flagged loudly rather than silently guessed:
- **Direct-entry J-tube geometry** (Fig. 13, no straight lead-in
  section). `full_pullin_analysis(direct_entry=True)` raises
  `NotImplementedError` -- the paper's one-paragraph description
  under-specifies this case, and two reconstructions were tried and
  both missed the paper's own stated result by more than the ~1-2%
  every other number here achieves. See
  `stage1_l1_over_d_direct_entry()`'s docstring for what was tried.
- **Fig. 6** (critical J-tube bend angle) -- given only as a chart,
  derived from a stability criterion not stated in closed form in the
  paper's body.
- **Multi-bend J-tubes** -- the paper explicitly defers this to "a
  forthcoming paper" -- see Ronalds (1992) below for the multi-bend
  extension this project does implement, though for a different
  (flexible-cable) physical regime.
- **Back-tension / pull-cable friction wrap**,
  `P_T = (P1+T)*exp(mu_c*psi) + W` -- a simple multiplicative step the
  caller can apply on top of `P1`/`P_III` from this module; not itself
  part of the J-tube mechanics being verified here.

### Verification approach

Every number is checked against the paper's own Appendix II worked
example (400 N/mm² yield riser, R/d=100, D/d=1.5, d=0.3 m, t=20 mm,
mu=0.3) -- see `results/appendix2_comparison.txt` for the full table.
Discrepancies sit at 0.0-1.7%, and are expected: **the paper's own
worked example reads two of its key quantities (alpha·Mbar/Mp,
beta·Mbar/Mp) off a hand-drawn chart (Fig. 12)** rather than evaluating
its own eqs. (7)-(8) directly. This module evaluates the equations,
which is more precise than a hand graph-read and is the only way to
compute results for parameter combinations the paper never plotted.

One genuine inconsistency was found in the source paper, not a
graph-reading artifact: the **Fig. 7 caption** (as transcribed from
the scan) defines `Po = pi*dm*t*sigma_o` using the mean diameter, but
the **Appendix II worked example's own arithmetic** only reproduces
the stated 7540 kN using the **outer** diameter `d` (mean diameter
gives 7037 kN, 6.7% off -- too large to be rounding). This module
follows the worked example, since it is independently verifiable and
a scanned figure-caption subscript is not. See
`full_pullin_analysis()`'s `Po` calculation for the full note.

Design-chart reproduction (Figs. 5, 7) matches the paper's plotted
shape and magnitude across the full R/d, D/d range, not just the one
worked-example point -- see `plots/`.

### GUI

`walker1983_gui.py` -- tkinter (stdlib) + an embedded matplotlib
design-curve plot. Opens pre-filled with the Appendix II parameters;
edit any field and press Compute (or Enter) to re-run. Shows a
paper-vs-module comparison column automatically, but only when the
current inputs match the Appendix II reference case exactly (there is
no paper-stated value to compare against for any other geometry). The
plot updates to show where the current case sits on the Stage I
pull-in-load-vs-R/d curve. Checking "Direct-entry geometry" surfaces
the `NotImplementedError` from above rather than silently computing
the standard case under a misleading label.

## Ronalds (1992) -- what's implemented

A different physical regime from Walker & Davies: the umbilical is
treated as flexible in bending, carrying load essentially as a cable
(friction-capstan tension growth, not elastic-plastic pipe bending).

| Piece | Status |
|---|---|
| Capstan/friction tension growth around one bend (eq. 1-3) | Verified against the textbook capstan equation and an independent ODE integration |
| `M = TL*y` moment identity (eq. 4) | Verified via direct statics for the simplest case |
| Multi-bend tension growth (eq. 10-11) | Verified: matches sequential single-bend application and the "total angle is what matters" invariant |
| Single-bend GUIDE support: eta, MB, max span moment (eq. 5-9) | Cross-checked against the paper's own eq. (6) numerical limit (`eta -> 2/pi`) |
| Single-bend FREE and FULL SUPPORT cases (Fig. 5 table) | FULL SUPPORT implemented exactly (paper states it unconditionally); FREE support's `MB` row only |
| Umbilical wrap-around-bend mechanics (eq. 12-18) | Implemented and sanity-checked |

**Not implemented**, flagged rather than guessed:
- The FREE-support-case `V`/`H` formulas and **all four deflection
  formulas** (`(EI)*Δx`, `(EI)*Δy`, both support cases) from the Fig. 5
  table -- the densest expressions in the source scan, with no
  independent numerical anchor to catch a transcription error against
  (unlike the guide-case `eta`, which eq. 6 does give a number for).
- Fig. 6's general multi-bend, multi-plane moment/reaction solve --
  the paper itself says this needs "a stress analysis program," not a
  closed form.
- Diaphragm pull-in (Fig. 8) and accidental pull-out (Fig. 9) -- reuse
  the same Fig. 5 machinery with different applied end forces, per the
  paper's own text; no new equations to implement.

### Verification approach -- weaker than Walker & Davies, and why

Ronalds (1992) gives **no worked numerical example anywhere in the
paper** -- everything is symbolic. This is a materially weaker basis
for verification than Walker & Davies, which pins every quantity
against a stated number. What this module relies on instead:

1. The one genuine numerical fact the paper does state: eq. (6)'s
   `eta -> 2/pi` limit, which the general guide-support formula is
   checked against and matches to 1e-12.
2. First-principles re-derivations that don't depend on trusting the
   OCR transcription at all -- the capstan equation checked against
   its standard textbook form and against direct ODE integration; the
   `M=TL*y` identity checked against plain statics on a bare
   quarter-circle.
3. Unconditional statements in the paper's prose (the full-support
   case) that need no formula transcription.
4. Physical monotonicity checks (friction only ever raises tension,
   etc).

Treat the guide-support `MB`/`eta` formulas as "transcribed and
internally consistent with the one number the paper gives," not
"independently verified" in the stronger sense the Walker & Davies
module achieves.

### GUI

`ronalds1992_gui.py` -- same tkinter + matplotlib convention as
`walker1983_gui.py`, with two tabs matching the module's two
independent calculation modes:

- **Multi-bend tension growth**: back tension, friction coefficient,
  and an editable list of bend angles in; tension after each bend and
  a tension-vs-cumulative-angle plot out. Opens pre-filled with the
  90°+8.5°+8.5° geometry used in Fan et al. (2013)'s multi-bend
  example, for comparison purposes only (Ronalds gives no worked
  example of its own -- see above).
- **Single-bend support response**: gamma/lambda/xi and a Free/Guide/
  Full-support selector in; eta, MB, and (for the Guide case) the max
  span moment out. Always shows the eq. (6) reference value (`2/pi`)
  alongside the current Guide-case result, since that is the one
  number this module's formula is actually checked against -- there is
  no equivalent "paper vs. module" comparison possible here the way
  `walker1983_gui.py` provides for its Appendix II case.

## Engineering note

`results/JTube_Method_Selection_Note.docx` -- a standalone, hand-off-
ready write-up of the "Are the two theories the same?" section above:
purpose, the two methods side by side, why cross-validation is
invalid, the wrap-vs-bridge test, the worked example with its
supporting figure, and a recommendation. Doesn't require reading this
README or the code to follow.

## License

MIT -- see [LICENSE](LICENSE). Free to use, modify, and redistribute,
with attribution.
