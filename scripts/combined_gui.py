# -*- coding: utf-8 -*-
"""
Combined desktop GUI: Walker & Davies (1983), Ronalds (1992), and the
cross-paper regime classifier, in one window with three tabs.

Each tab is functionally equivalent to its standalone counterpart
(walker1983_gui.py, ronalds1992_gui.py) -- this file does not change
either module's calculation logic, only re-hosts the same input/
compute/display pattern inside shared tabs so a single geometry can
flow between them without retyping numbers. The standalone GUIs remain
untouched and fully usable on their own.

The one genuinely NEW piece of behaviour: the Regime Classifier tab
has a "Use Walker tab's P1" button that pulls the currently-computed
Stage I pull load straight from the Walker tab into the classifier's
pull-tension field -- the same manual step used in the terminal demo
that produced the ~11.5x margin finding documented in README.md.

Run:
    python combined_gui.py
"""

import math
import os
import sys
import tkinter as tk
from tkinter import ttk

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import walker1983_jtube as jt
import ronalds1992_jtube as rn
import regime_classifier as rc


WALKER_DEFAULTS = dict(
    sigma_o_MPa='400', E_GPa='200', R_m='30', D_m='0.45',
    d_m='0.30', t_mm='20', mu='0.30',
)
WALKER_REFERENCE = {
    'Po (kN)': 7540.0, 'Mbar/Mp': 0.97, 'Stage I: L/d': 32.7,
    'Stage I: l1/d': 12.6, 'Stage I: tau (deg)': 7.1,
    'Stage I: P1/Po': 2.2e-2, 'Stage I: P1 (kN)': 166.0,
    'Stage II: l2/d': 34.5, 'Stage II: P2/Po': 1.27e-2,
    'Stage III: Rp/d': 166.7, "Stage III: P'/Po (peak)": 0.92e-2,
    'Stage III: P3/Po': 2.19e-2, 'Stage III: P3 (kN)': None,
}
DEFAULT_BENDS = [90.0, 8.5, 8.5]


# =========================================================================
# Tab 1: Walker & Davies (1983)
# =========================================================================

class WalkerTab(ttk.Frame):

    def __init__(self, parent):
        super().__init__(parent)
        self.last_P1_kN = None   # exposed for the Regime Classifier tab
        self._build()
        self.compute()

    def _build(self):
        self.columnconfigure(1, weight=1)
        self.rowconfigure(1, weight=1)

        inputs = ttk.LabelFrame(self, text='Inputs  (defaults = Appendix II worked example)')
        inputs.grid(row=0, column=0, sticky='nsew', padx=8, pady=8)

        fields = [
            ('sigma_o_MPa', 'Riser yield stress, sigma_o (MPa)'),
            ('E_GPa', "Young's modulus, E (GPa)"),
            ('R_m', 'J-tube bend radius, R (m)'),
            ('D_m', 'J-tube inside diameter, D (m)'),
            ('d_m', 'Riser outside diameter, d (m)'),
            ('t_mm', 'Riser wall thickness, t (mm)'),
            ('mu', 'Coefficient of friction, mu'),
        ]
        self.vars = {}
        for i, (key, label) in enumerate(fields):
            ttk.Label(inputs, text=label).grid(row=i, column=0, sticky='w', padx=4, pady=3)
            var = tk.StringVar(value=WALKER_DEFAULTS[key])
            entry = ttk.Entry(inputs, textvariable=var, width=12)
            entry.grid(row=i, column=1, padx=4, pady=3)
            entry.bind('<Return>', lambda e: self.compute())
            self.vars[key] = var

        self.direct_entry_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(inputs, text='Direct-entry geometry (Fig. 13)\n'
                                     '[not implemented -- see below]',
                        variable=self.direct_entry_var).grid(
            row=len(fields), column=0, columnspan=2, sticky='w', padx=4, pady=(8, 3))

        ttk.Button(inputs, text='Compute', command=self.compute).grid(
            row=len(fields) + 1, column=0, columnspan=2, pady=(10, 4), sticky='ew')
        ttk.Button(inputs, text='Reset to Appendix II defaults',
                  command=self.reset_defaults).grid(
            row=len(fields) + 2, column=0, columnspan=2, pady=2, sticky='ew')

        self.status_var = tk.StringVar(value='')
        ttk.Label(inputs, textvariable=self.status_var, foreground='#97281f',
                 wraplength=220, justify='left').grid(
            row=len(fields) + 3, column=0, columnspan=2, sticky='w', padx=4, pady=(8, 4))

        results = ttk.LabelFrame(self, text='Results  (this module vs. the paper-stated value, where available)')
        results.grid(row=0, column=1, sticky='nsew', padx=8, pady=8)

        columns = ('quantity', 'value', 'paper', 'diff')
        self.tree = ttk.Treeview(results, columns=columns, show='headings', height=14)
        headers = {'quantity': 'Quantity', 'value': 'This module',
                  'paper': 'Paper (Appx. II)', 'diff': '% diff'}
        widths = {'quantity': 200, 'value': 100, 'paper': 110, 'diff': 70}
        for c in columns:
            self.tree.heading(c, text=headers[c])
            self.tree.column(c, width=widths[c], anchor='center' if c != 'quantity' else 'w')
        self.tree.pack(fill='both', expand=True, padx=4, pady=4)

        plot_frame = ttk.LabelFrame(self, text='Stage I pull-in load vs. R/d')
        plot_frame.grid(row=1, column=0, columnspan=2, sticky='nsew', padx=8, pady=(0, 8))
        self.fig, self.ax = plt.subplots(figsize=(9, 2.8))
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.get_tk_widget().pack(fill='both', expand=True, padx=4, pady=4)

    def reset_defaults(self):
        for key, val in WALKER_DEFAULTS.items():
            self.vars[key].set(val)
        self.direct_entry_var.set(False)
        self.compute()

    def _read_inputs(self):
        vals = {}
        labels = {
            'sigma_o_MPa': 'yield stress', 'E_GPa': "Young's modulus",
            'R_m': 'bend radius', 'D_m': 'J-tube ID', 'd_m': 'riser OD',
            't_mm': 'wall thickness', 'mu': 'friction coefficient',
        }
        for key, var in self.vars.items():
            raw = var.get().strip()
            try:
                vals[key] = float(raw)
            except ValueError:
                raise ValueError(f'"{raw}" is not a valid number for {labels[key]}.')
        for key in ('sigma_o_MPa', 'E_GPa', 'R_m', 'D_m', 'd_m', 't_mm'):
            if vals[key] <= 0:
                raise ValueError(f'{labels[key]} must be positive (got {vals[key]}).')
        if not (0.0 <= vals['mu'] < 2.0):
            raise ValueError(f"friction coefficient mu={vals['mu']} is outside a "
                             f"physically sensible range [0, 2).")
        if vals['d_m'] >= vals['D_m']:
            raise ValueError(f"riser OD (d={vals['d_m']} m) must be less than "
                             f"J-tube ID (D={vals['D_m']} m).")
        if vals['t_mm'] / 1000.0 >= vals['d_m'] / 2.0:
            raise ValueError(f"wall thickness (t={vals['t_mm']} mm) is not physically "
                             f"possible for riser OD d={vals['d_m']} m.")
        return vals

    def compute(self):
        self.status_var.set('')
        try:
            vals = self._read_inputs()
        except ValueError as e:
            self.status_var.set(str(e))
            self._clear_results()
            self.last_P1_kN = None
            return

        sigma_o = vals['sigma_o_MPa'] * 1e6
        E = vals['E_GPa'] * 1e9
        R, D, d, t, mu = vals['R_m'], vals['D_m'], vals['d_m'], vals['t_mm'] / 1000.0, vals['mu']
        direct_entry = self.direct_entry_var.get()

        try:
            result = jt.full_pullin_analysis(R, D, d, t, sigma_o, E, mu, direct_entry=direct_entry)
        except NotImplementedError:
            self.status_var.set(
                'Direct-entry geometry is not implemented (see walker1983_jtube.py '
                'docstring). Showing the standard lead-in-section case instead.'
            )
            try:
                result = jt.full_pullin_analysis(R, D, d, t, sigma_o, E, mu, direct_entry=False)
            except ValueError as e2:
                self._clear_results()
                self.status_var.set(f'Non-physical inputs: {e2}')
                self.last_P1_kN = None
                return
        except ValueError as e:
            self._clear_results()
            self.status_var.set(f'Non-physical result for these inputs: {e}')
            self.last_P1_kN = None
            return

        self._populate_results(result, vals)
        self._update_plot(result, vals)
        self.last_P1_kN = result['P1'] / 1e3

    def _clear_results(self):
        for row in self.tree.get_children():
            self.tree.delete(row)

    def _populate_results(self, result, vals):
        self._clear_results()
        s1, s2, s3 = result['stage1'], result['stage2'], result['stage3']
        rows = [
            ('Po (kN)', result['Po'] / 1e3), ('Mbar/Mp', s1['Mbar_Mp']),
            ('Stage I: L/d', s1['L_over_d']), ('Stage I: l1/d', s1['l1_over_d']),
            ('Stage I: tau (deg)', s1['tau_deg']), ('Stage I: P1/Po', s1['P1_over_Po']),
            ('Stage I: P1 (kN)', result['P1'] / 1e3), ('Stage II: l2/d', s2['l2_over_d']),
            ('Stage II: P2/Po', s2['P2_over_Po']), ('Stage III: Rp/d', s3['Rp_over_d']),
            ("Stage III: P'/Po (peak)", s3['P_prime_over_Po']),
            ('Stage III: P3/Po', s3['P3_over_Po']), ('Stage III: P3 (kN)', result['P3'] / 1e3),
        ]
        is_reference_case = (
            abs(vals['sigma_o_MPa'] - 400) < 1e-6 and abs(vals['E_GPa'] - 200) < 1e-6
            and abs(vals['R_m'] - 30) < 1e-6 and abs(vals['D_m'] - 0.45) < 1e-6
            and abs(vals['d_m'] - 0.30) < 1e-6 and abs(vals['t_mm'] - 20) < 1e-6
            and abs(vals['mu'] - 0.30) < 1e-6
        )
        for name, value in rows:
            val_str = f'{value:.4g}'
            paper_val = WALKER_REFERENCE.get(name)
            if is_reference_case and paper_val is not None:
                paper_str = f'{paper_val:.4g}'
                diff_str = f'{100.0*(value-paper_val)/paper_val:+.1f}%'
            else:
                paper_str = '' if not is_reference_case else '-- (not stated)'
                diff_str = ''
            self.tree.insert('', 'end', values=(name, val_str, paper_str, diff_str))

    def _update_plot(self, result, vals):
        self.ax.clear()
        Dd = vals['D_m'] / vals['d_m']
        sigma_o_E = vals['sigma_o_MPa'] * 1e6 / (vals['E_GPa'] * 1e9)
        Rd_range = np.linspace(20.0, 300.0, 141)
        Rds, P1s = [], []
        for Rd in Rd_range:
            try:
                r = jt.stage1_pullin(Rd, Dd, sigma_o=sigma_o_E, E=1.0, mu=vals['mu'])
                Rds.append(Rd)
                P1s.append(r['P1_over_Po'] * 1e3)
            except ValueError:
                continue
        self.ax.plot(Rds, P1s, color='#1a5c78', label=f'D/d={Dd:.3g}, mu={vals["mu"]:.3g}')
        Rd_current = vals['R_m'] / vals['d_m']
        P1_current = result['stage1']['P1_over_Po'] * 1e3
        self.ax.plot([Rd_current], [P1_current], 'o', color='#97281f', markersize=9, zorder=5,
                     label=f'Current: R/d={Rd_current:.1f}, P1/Po={P1_current:.2f}e-3')
        self.ax.set_xlabel('R/d')
        self.ax.set_ylabel('(P1/Po) x 10^3')
        self.ax.legend(fontsize=8, loc='upper right')
        self.ax.grid(True, alpha=0.3)
        self.fig.tight_layout()
        self.canvas.draw()


# =========================================================================
# Tab 2: Ronalds (1992)
# =========================================================================

class RonaldsTab(ttk.Frame):

    def __init__(self, parent):
        super().__init__(parent)
        self._build()
        self.compute_multibend()
        self.compute_support()

    def _build(self):
        inner = ttk.Notebook(self)
        inner.pack(fill='both', expand=True)
        self.tab_multibend = ttk.Frame(inner)
        self.tab_support = ttk.Frame(inner)
        inner.add(self.tab_multibend, text='Multi-bend tension growth')
        inner.add(self.tab_support, text='Single-bend support response')
        self._build_multibend_tab()
        self._build_support_tab()

    def _build_multibend_tab(self):
        tab = self.tab_multibend
        tab.columnconfigure(1, weight=1)
        tab.rowconfigure(1, weight=1)

        inputs = ttk.LabelFrame(tab, text='Inputs')
        inputs.grid(row=0, column=0, sticky='nsew', padx=8, pady=8)

        ttk.Label(inputs, text='Back tension T0 (kN)').grid(row=0, column=0, sticky='w', padx=4, pady=3)
        self.T0_var = tk.StringVar(value='100')
        ttk.Entry(inputs, textvariable=self.T0_var, width=10).grid(row=0, column=1, padx=4, pady=3)

        ttk.Label(inputs, text='Friction coefficient, mu').grid(row=1, column=0, sticky='w', padx=4, pady=3)
        self.mu_mb_var = tk.StringVar(value='0.30')
        ttk.Entry(inputs, textvariable=self.mu_mb_var, width=10).grid(row=1, column=1, padx=4, pady=3)

        ttk.Label(inputs, text='Bend angles (deg), one per line:').grid(
            row=2, column=0, columnspan=2, sticky='w', padx=4, pady=(10, 2))
        self.bends_text = tk.Text(inputs, width=14, height=6)
        self.bends_text.grid(row=3, column=0, columnspan=2, padx=4, pady=2)
        self.bends_text.insert('1.0', '\n'.join(str(a) for a in DEFAULT_BENDS))

        ttk.Button(inputs, text='Compute', command=self.compute_multibend).grid(
            row=4, column=0, columnspan=2, pady=(10, 4), sticky='ew')

        self.mb_status_var = tk.StringVar(value='')
        ttk.Label(inputs, textvariable=self.mb_status_var, foreground='#97281f',
                 wraplength=180, justify='left').grid(
            row=5, column=0, columnspan=2, sticky='w', padx=4, pady=(8, 4))

        results = ttk.LabelFrame(tab, text='Tension after each bend')
        results.grid(row=0, column=1, sticky='nsew', padx=8, pady=8)
        columns = ('bend', 'angle_deg', 'tension_kN', 'ratio')
        self.mb_tree = ttk.Treeview(results, columns=columns, show='headings', height=8)
        headers = {'bend': 'Bend #', 'angle_deg': 'Angle (deg)',
                  'tension_kN': 'Tension (kN)', 'ratio': 'T / T0'}
        for c in columns:
            self.mb_tree.heading(c, text=headers[c])
            self.mb_tree.column(c, width=100, anchor='center')
        self.mb_tree.pack(fill='both', expand=True, padx=4, pady=4)

        plot_frame = ttk.LabelFrame(tab, text='Tension vs. cumulative bend angle')
        plot_frame.grid(row=1, column=0, columnspan=2, sticky='nsew', padx=8, pady=(0, 8))
        tab.rowconfigure(1, weight=1)
        self.mb_fig, self.mb_ax = plt.subplots(figsize=(8, 2.2))
        self.mb_canvas = FigureCanvasTkAgg(self.mb_fig, master=plot_frame)
        self.mb_canvas.get_tk_widget().pack(fill='both', expand=True, padx=4, pady=4)

    def compute_multibend(self):
        self.mb_status_var.set('')
        raw_T0 = self.T0_var.get().strip()
        try:
            T0 = float(raw_T0) * 1e3
        except ValueError:
            self.mb_status_var.set(f'"{raw_T0}" is not a valid number for back tension T0.')
            return
        raw_mu = self.mu_mb_var.get().strip()
        try:
            mu = float(raw_mu)
        except ValueError:
            self.mb_status_var.set(f'"{raw_mu}" is not a valid number for mu.')
            return
        try:
            if T0 <= 0:
                raise ValueError('back tension T0 must be positive.')
            if not (0.0 <= mu < 2.0):
                raise ValueError(f'mu={mu} is outside a physically sensible range [0, 2).')
            raw_lines = [l.strip() for l in self.bends_text.get('1.0', 'end').splitlines() if l.strip()]
            angles_deg = []
            for line in raw_lines:
                try:
                    angles_deg.append(float(line))
                except ValueError:
                    raise ValueError(f'"{line}" is not a valid bend angle.')
            if not angles_deg:
                raise ValueError('enter at least one bend angle.')
            for a in angles_deg:
                if not (0.0 <= a <= 180.0):
                    raise ValueError(f'bend angle {a} deg is outside [0, 180].')
        except ValueError as e:
            self.mb_status_var.set(str(e))
            return

        alphas = [math.radians(a) for a in angles_deg]
        tensions = rn.tension_after_multiple_bends(T0, mu, alphas)
        for row in self.mb_tree.get_children():
            self.mb_tree.delete(row)
        for i, (a_deg, T) in enumerate(zip(angles_deg, tensions[1:]), start=1):
            self.mb_tree.insert('', 'end', values=(i, f'{a_deg:.2f}', f'{T/1e3:.3f}', f'{T/T0:.4f}'))

        self.mb_ax.clear()
        cum_angle = np.concatenate([[0.0], np.cumsum(angles_deg)])
        self.mb_ax.plot(cum_angle, np.array(tensions) / 1e3, marker='o', color='#1a5c78')
        self.mb_ax.set_xlabel('Cumulative bend angle (deg)')
        self.mb_ax.set_ylabel('Tension (kN)')
        self.mb_ax.grid(True, alpha=0.3)
        self.mb_fig.tight_layout()
        self.mb_canvas.draw()

    def _build_support_tab(self):
        tab = self.tab_support
        tab.columnconfigure(1, weight=1)
        tab.rowconfigure(1, weight=1)

        inputs = ttk.LabelFrame(tab, text='Inputs')
        inputs.grid(row=0, column=0, sticky='nsew', padx=8, pady=8)

        ttk.Label(inputs, text='Support type at A').grid(row=0, column=0, sticky='w', padx=4, pady=3)
        self.support_var = tk.StringVar(value='Guide')
        cb = ttk.Combobox(inputs, textvariable=self.support_var, state='readonly',
                          values=['Free', 'Guide', 'Full support'], width=12)
        cb.grid(row=0, column=1, padx=4, pady=3)
        self.support_var.trace_add('write', lambda *a: self.compute_support())

        ttk.Label(inputs, text='gamma (below bend, /r)').grid(row=1, column=0, sticky='w', padx=4, pady=3)
        self.gamma_var = tk.StringVar(value='0.0')
        ttk.Entry(inputs, textvariable=self.gamma_var, width=10).grid(row=1, column=1, padx=4, pady=3)

        ttk.Label(inputs, text='lambda (B to D, /r)').grid(row=2, column=0, sticky='w', padx=4, pady=3)
        self.lam_var = tk.StringVar(value='0.0')
        ttk.Entry(inputs, textvariable=self.lam_var, width=10).grid(row=2, column=1, padx=4, pady=3)

        ttk.Label(inputs, text='xi (bend top to B, /r)').grid(row=3, column=0, sticky='w', padx=4, pady=3)
        self.xi_var = tk.StringVar(value='0.0')
        ttk.Entry(inputs, textvariable=self.xi_var, width=10).grid(row=3, column=1, padx=4, pady=3)

        ttk.Button(inputs, text='Compute', command=self.compute_support).grid(
            row=4, column=0, columnspan=2, pady=(10, 4), sticky='ew')

        self.sup_status_var = tk.StringVar(value='')
        ttk.Label(inputs, textvariable=self.sup_status_var, foreground='#97281f',
                 wraplength=180, justify='left').grid(
            row=5, column=0, columnspan=2, sticky='w', padx=4, pady=(8, 4))

        results = ttk.LabelFrame(tab, text='Result  (normalised by TL and r, where applicable)')
        results.grid(row=0, column=1, sticky='nsew', padx=8, pady=8)
        columns = ('quantity', 'value', 'reference')
        self.sup_tree = ttk.Treeview(results, columns=columns, show='headings', height=8)
        headers = {'quantity': 'Quantity', 'value': 'Value', 'reference': "Paper's eq. (6) reference"}
        widths = {'quantity': 180, 'value': 100, 'reference': 180}
        for c in columns:
            self.sup_tree.heading(c, text=headers[c])
            self.sup_tree.column(c, width=widths[c], anchor='w' if c != 'value' else 'center')
        self.sup_tree.pack(fill='both', expand=True, padx=4, pady=4)

    def compute_support(self):
        self.sup_status_var.set('')
        fields = [('gamma', self.gamma_var), ('lambda', self.lam_var), ('xi', self.xi_var)]
        values = {}
        for name, var in fields:
            raw = var.get().strip()
            try:
                values[name] = float(raw)
            except ValueError:
                self.sup_status_var.set(f'"{raw}" is not a valid number for {name}.')
                return
        gamma, lam, xi = values['gamma'], values['lambda'], values['xi']
        if gamma < 0 or lam < 0 or xi < 0:
            self.sup_status_var.set("gamma, lambda, xi must all be >= 0 (the paper's own stated restriction).")
            return

        support = self.support_var.get()
        for row in self.sup_tree.get_children():
            self.sup_tree.delete(row)

        if support == 'Full support':
            r = rn.full_support_response()
            rows = [('eta = V/TL', r['eta'], ''), ('MB/(TL*r)', r['MB_over_TLr'], ''),
                    ('H/TL', r['H_over_TL'], '')]
        elif support == 'Free':
            MB = rn.free_support_MB_over_TLr(xi)
            rows = [('MB/(TL*r)', MB, ''), ('eta = V/TL', 0.0, '(V=0 for the free case)'),
                    ('H/TL', 0.0, '(H=0 for the free case)')]
        else:
            r = rn.guide_support_response(gamma, lam, xi)
            ref_note = f"matches eq.(6): 2/pi = {2/math.pi:.6f}" if (gamma == 0 and lam == 0 and xi == 0) else ''
            rows = [('eta = V/TL', r['eta'], ref_note), ('MB/(TL*r)', r['MB_over_TLr'], ''),
                    ('theta at max Ms (deg)', r['theta_max_deg'], ''), ('MS/(TL*r)', r['MS_over_TLr'], '')]

        for name, value, note in rows:
            self.sup_tree.insert('', 'end', values=(name, f'{value:.5f}', note))


# =========================================================================
# Tab 3: Regime Classifier
# =========================================================================

class RegimeTab(ttk.Frame):

    def __init__(self, parent, walker_tab):
        super().__init__(parent)
        self.walker_tab = walker_tab   # for "Use Walker tab's P1"
        self._build()
        self.compute()

    def _build(self):
        self.columnconfigure(1, weight=1)
        self.rowconfigure(1, weight=1)

        inputs = ttk.LabelFrame(self, text='Inputs')
        inputs.grid(row=0, column=0, sticky='nsew', padx=8, pady=8)

        fields = [
            ('OD_m', 'Outer diameter, OD (m)', '0.30'),
            ('t_mm', 'Wall thickness, t (mm)', '20'),
            ('E_GPa', "Young's modulus, E (GPa)", '200'),
            ('r_m', 'J-tube bend radius, r (m)', '30'),
            ('D_m', 'J-tube inside diameter, D (m)', '0.45'),
        ]
        self.vars = {}
        for i, (key, label, default) in enumerate(fields):
            ttk.Label(inputs, text=label).grid(row=i, column=0, sticky='w', padx=4, pady=3)
            var = tk.StringVar(value=default)
            ttk.Entry(inputs, textvariable=var, width=12).grid(row=i, column=1, padx=4, pady=3)
            self.vars[key] = var

        n = len(fields)
        ttk.Label(inputs, text='Estimated pull tension (kN)\n(leave blank to see the '
                              'threshold only)').grid(
            row=n, column=0, columnspan=2, sticky='w', padx=4, pady=(10, 2))
        self.T_var = tk.StringVar(value='167.2')
        ttk.Entry(inputs, textvariable=self.T_var, width=12).grid(row=n + 1, column=0, padx=4, pady=2, sticky='w')

        ttk.Button(inputs, text='Use Walker tab\'s P1',
                  command=self._use_walker_P1).grid(row=n + 2, column=0, columnspan=2, pady=(4, 2), sticky='ew')

        ttk.Button(inputs, text='Compute', command=self.compute).grid(
            row=n + 3, column=0, columnspan=2, pady=(10, 4), sticky='ew')
        ttk.Button(inputs, text='Load Walker Appendix II riser',
                  command=self._load_walker_example).grid(row=n + 4, column=0, columnspan=2, pady=2, sticky='ew')
        ttk.Button(inputs, text='Load placeholder flexible umbilical',
                  command=self._load_flexible_example).grid(row=n + 5, column=0, columnspan=2, pady=2, sticky='ew')

        self.status_var = tk.StringVar(value='')
        ttk.Label(inputs, textvariable=self.status_var, foreground='#97281f',
                 wraplength=220, justify='left').grid(
            row=n + 6, column=0, columnspan=2, sticky='w', padx=4, pady=(8, 4))

        results = ttk.LabelFrame(self, text='Classification')
        results.grid(row=0, column=1, sticky='nsew', padx=8, pady=8)

        self.regime_var = tk.StringVar(value='')
        self.regime_label = ttk.Label(results, textvariable=self.regime_var,
                                      font=('TkDefaultFont', 13, 'bold'), wraplength=380)
        self.regime_label.pack(anchor='w', padx=10, pady=(10, 6))

        columns = ('quantity', 'value')
        self.tree = ttk.Treeview(results, columns=columns, show='headings', height=6)
        self.tree.heading('quantity', text='Quantity')
        self.tree.heading('value', text='Value')
        self.tree.column('quantity', width=220, anchor='w')
        self.tree.column('value', width=140, anchor='center')
        self.tree.pack(fill='both', expand=True, padx=10, pady=6)

        note = ('The wrap-vs-bridge threshold is Ronalds (1992) eq. (15): the '
               'tension above which this object actually conforms continuously '
               'to the bend, rather than touching it at discrete points as '
               "Walker & Davies' method assumes. See README.md \"Are the two "
               'theories the same?" for what this does and does not prove.')
        ttk.Label(results, text=note, foreground='#697680', wraplength=380,
                 justify='left', font=('TkDefaultFont', 8)).pack(anchor='w', padx=10, pady=(6, 10))

        plot_frame = ttk.LabelFrame(self, text='Wrap threshold vs. bend radius r  (at the current OD, t, E, D)')
        plot_frame.grid(row=1, column=0, columnspan=2, sticky='nsew', padx=8, pady=(0, 8))
        self.fig, self.ax = plt.subplots(figsize=(9, 2.8))
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.get_tk_widget().pack(fill='both', expand=True, padx=4, pady=4)

    def _use_walker_P1(self):
        if self.walker_tab.last_P1_kN is None:
            self.status_var.set("Walker tab has no valid computed P1 yet -- "
                                "switch to it, fix any input errors, and try again.")
            return
        self.T_var.set(f'{self.walker_tab.last_P1_kN:.4g}')
        self.compute()

    def _load_walker_example(self):
        self.vars['OD_m'].set('0.30')
        self.vars['t_mm'].set('20')
        self.vars['E_GPa'].set('200')
        self.vars['r_m'].set('30')
        self.vars['D_m'].set('0.45')
        self.T_var.set('167.2')
        self.compute()

    def _load_flexible_example(self):
        self.vars['OD_m'].set('0.10')
        self.vars['t_mm'].set('10')
        self.vars['E_GPa'].set('1')
        self.vars['r_m'].set('12')
        self.vars['D_m'].set('0.20')
        self.T_var.set('50')
        self.compute()

    def compute(self):
        self.status_var.set('')
        labels = {'OD_m': 'outer diameter', 't_mm': 'wall thickness',
                 'E_GPa': "Young's modulus", 'r_m': 'bend radius', 'D_m': 'J-tube ID'}
        vals = {}
        for key, var in self.vars.items():
            raw = var.get().strip()
            try:
                vals[key] = float(raw)
            except ValueError:
                self.status_var.set(f'"{raw}" is not a valid number for {labels[key]}.')
                return
            if vals[key] <= 0:
                self.status_var.set(f'{labels[key]} must be positive (got {vals[key]}).')
                return

        raw_T = self.T_var.get().strip()
        T = None
        if raw_T:
            try:
                T = float(raw_T) * 1e3
            except ValueError:
                self.status_var.set(f'"{raw_T}" is not a valid number for pull tension.')
                return
            if T <= 0:
                self.status_var.set(f'pull tension must be positive (got {raw_T}).')
                return

        OD, t, E, r, D = vals['OD_m'], vals['t_mm'] / 1000.0, vals['E_GPa'] * 1e9, vals['r_m'], vals['D_m']
        if OD >= D:
            self.status_var.set(f'outer diameter (OD={OD} m) must be less than J-tube ID (D={D} m).')
            return
        if t >= OD / 2.0:
            self.status_var.set(f'wall thickness (t={t*1000:.1f} mm) is not physically possible for OD={OD} m.')
            return

        try:
            result = rc.classify_pullin_regime(OD, t, E, r, D, pullin_tension=T)
        except ValueError as e:
            self.status_var.set(str(e))
            return

        colors = {
            'rigid (Walker-type): discrete contact, nowhere near wrapping': '#1a5c78',
            'flexible (Ronalds-type): well past the wrap threshold': '#15624b',
            'marginal -- within 2x of the wrap threshold, check both': '#8f5307',
            'unknown (no pull tension given)': '#697680',
        }
        self.regime_var.set(result['regime'])
        self.regime_label.configure(foreground=colors.get(result['regime'], '#151b20'))

        for row in self.tree.get_children():
            self.tree.delete(row)
        self.tree.insert('', 'end', values=('EI', f"{result['EI']/1e6:.2f} MN*m^2"))
        self.tree.insert('', 'end', values=('Wrap threshold tension', f"{result['threshold_T']/1e3:.2f} kN"))
        if result['pullin_tension'] is not None:
            self.tree.insert('', 'end', values=('Pull tension (input)', f"{result['pullin_tension']/1e3:.2f} kN"))
            self.tree.insert('', 'end', values=('Margin (threshold / pull)', f"{result['margin']:.2f}x"))

        self._update_plot(OD, t, E, D, T)

    def _update_plot(self, OD, t, E, D, T):
        self.ax.clear()
        r_range = np.linspace(max(2.0, D * 1.5), 100.0, 100)
        thresholds = []
        for r in r_range:
            try:
                thresholds.append(rc.wrap_threshold_tension(OD, t, E, r, D) / 1e3)
            except ValueError:
                thresholds.append(np.nan)
        self.ax.plot(r_range, thresholds, color='#1a5c78', label='Wrap threshold (kN)')
        if T is not None:
            self.ax.axhline(T / 1e3, color='#97281f', linestyle=':', label=f'Current pull tension: {T/1e3:.1f} kN')
        self.ax.set_xlabel('J-tube bend radius r (m)')
        self.ax.set_ylabel('Threshold tension (kN)')
        self.ax.set_yscale('log')
        self.ax.legend(fontsize=8)
        self.ax.grid(True, alpha=0.3, which='both')
        self.fig.tight_layout()
        self.canvas.draw()


# =========================================================================
# Main window
# =========================================================================

class CombinedGUI(tk.Tk):

    def __init__(self):
        super().__init__()
        self.title('J-Tube Pull-In Analysis -- Walker & Davies / Ronalds / Regime Classifier')
        self.geometry('1260x760')
        self.minsize(1040, 640)
        self.protocol("WM_DELETE_WINDOW", self._on_close)

        notebook = ttk.Notebook(self)
        notebook.pack(fill='both', expand=True, padx=6, pady=6)

        walker_tab = WalkerTab(notebook)
        ronalds_tab = RonaldsTab(notebook)
        regime_tab = RegimeTab(notebook, walker_tab=walker_tab)

        notebook.add(walker_tab, text='Walker & Davies (1983)')
        notebook.add(ronalds_tab, text='Ronalds (1992)')
        notebook.add(regime_tab, text='Regime Classifier')

    def _on_close(self):
        """
        Explicitly close every embedded matplotlib figure and force the
        process to exit. Without this, closing the window can leave the
        Python process running in the background -- pyplot keeps its
        own global reference to every figure created via plt.subplots(),
        and this GUI creates four of them (one per tab, two on the
        Ronalds tab), so plain self.destroy() is not reliably enough to
        let the interpreter shut down on its own.
        """
        plt.close('all')
        self.destroy()
        os._exit(0)


def main():
    app = CombinedGUI()
    app.mainloop()


if __name__ == '__main__':
    main()
