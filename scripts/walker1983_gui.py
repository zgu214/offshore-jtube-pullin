# -*- coding: utf-8 -*-
"""
Desktop GUI for walker1983_jtube.py.

Uses tkinter (Python standard library -- no extra install needed,
matching this project's "Python 3, stdlib only" convention for the
core library) plus matplotlib for the embedded design-curve plot,
which reproduce_figures.py and the test suite already depend on.

Opens pre-filled with the Appendix II worked-example parameters (see
walker1983_jtube.py's module docstring), so the first thing you see on
launch is the reference case this module was verified against. Edit
any field and press "Compute" (or Enter) to re-run the analysis.

Run:
    python walker1983_gui.py
"""

import os
import sys
import tkinter as tk
from tkinter import ttk, messagebox

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import walker1983_jtube as jt


# Appendix II worked-example values -- the reference case this module
# was verified against (see walker1983_jtube.py docstring / README).
DEFAULTS = dict(
    sigma_o_MPa='400',
    E_GPa='200',
    R_m='30',
    D_m='0.45',
    d_m='0.30',
    t_mm='20',
    mu='0.30',
)

# Paper-stated reference values for the comparison column (results/appendix2_comparison.txt)
REFERENCE = {
    'Po (kN)': 7540.0,
    'Mbar/Mp': 0.97,
    'Stage I: L/d': 32.7,
    'Stage I: l1/d': 12.6,
    'Stage I: tau (deg)': 7.1,
    'Stage I: P1/Po': 2.2e-2,
    'Stage I: P1 (kN)': 166.0,
    'Stage II: l2/d': 34.5,
    'Stage II: P2/Po': 1.27e-2,
    'Stage III: Rp/d': 166.7,
    "Stage III: P'/Po (peak)": 0.92e-2,
    'Stage III: P3/Po': 2.19e-2,
    'Stage III: P3 (kN)': None,   # not separately stated by the paper
}


class WalkerGUI(tk.Tk):

    def __init__(self):
        super().__init__()
        self.title('Walker & Davies (1983) -- J-Tube Pull-In Calculator')
        self.geometry('1180x680')
        self.minsize(980, 600)
        self.protocol("WM_DELETE_WINDOW", self._on_close)

        self._build_input_panel()
        self._build_results_panel()
        self._build_plot_panel()

        self.compute()  # populate with the Appendix II defaults on launch

    def _on_close(self):
        """
        Explicitly close the embedded matplotlib figure and force the
        process to exit. Without this, closing the window can leave the
        Python process running in the background -- pyplot keeps its own
        global reference to the figure created via plt.subplots(), and
        plain self.destroy() is not reliably enough to let the
        interpreter shut down on its own.
        """
        plt.close('all')
        self.destroy()
        os._exit(0)

    # -----------------------------------------------------------
    # Layout
    # -----------------------------------------------------------

    def _build_input_panel(self):
        frame = ttk.LabelFrame(self, text='Inputs  (defaults = Appendix II worked example)')
        frame.grid(row=0, column=0, sticky='nsew', padx=8, pady=8)

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
            ttk.Label(frame, text=label).grid(row=i, column=0, sticky='w', padx=4, pady=3)
            var = tk.StringVar(value=DEFAULTS[key])
            entry = ttk.Entry(frame, textvariable=var, width=12)
            entry.grid(row=i, column=1, padx=4, pady=3)
            entry.bind('<Return>', lambda e: self.compute())
            self.vars[key] = var

        self.direct_entry_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(frame, text='Direct-entry geometry (Fig. 13)\n'
                                     '[not implemented -- see below]',
                        variable=self.direct_entry_var).grid(
            row=len(fields), column=0, columnspan=2, sticky='w', padx=4, pady=(8, 3))

        ttk.Button(frame, text='Compute', command=self.compute).grid(
            row=len(fields) + 1, column=0, columnspan=2, pady=(10, 4), sticky='ew')
        ttk.Button(frame, text='Reset to Appendix II defaults',
                  command=self.reset_defaults).grid(
            row=len(fields) + 2, column=0, columnspan=2, pady=2, sticky='ew')

        self.status_var = tk.StringVar(value='')
        ttk.Label(frame, textvariable=self.status_var, foreground='#97281f',
                 wraplength=220, justify='left').grid(
            row=len(fields) + 3, column=0, columnspan=2, sticky='w', padx=4, pady=(8, 4))

        self.rowconfigure(0, weight=0)
        self.columnconfigure(1, weight=1)

    def _build_results_panel(self):
        frame = ttk.LabelFrame(self, text='Results  (this module vs. the paper-stated value, where available)')
        frame.grid(row=0, column=1, sticky='nsew', padx=8, pady=8)
        self.columnconfigure(1, weight=2)

        columns = ('quantity', 'value', 'paper', 'diff')
        self.tree = ttk.Treeview(frame, columns=columns, show='headings', height=16)
        headers = {'quantity': 'Quantity', 'value': 'This module',
                  'paper': 'Paper (Appx. II)', 'diff': '% diff'}
        widths = {'quantity': 220, 'value': 110, 'paper': 120, 'diff': 80}
        for c in columns:
            self.tree.heading(c, text=headers[c])
            self.tree.column(c, width=widths[c], anchor='center' if c != 'quantity' else 'w')
        self.tree.pack(fill='both', expand=True, padx=4, pady=4)

    def _build_plot_panel(self):
        frame = ttk.LabelFrame(self, text='Stage I pull-in load vs. R/d  (D/d fixed at the current input)')
        frame.grid(row=1, column=0, columnspan=2, sticky='nsew', padx=8, pady=(0, 8))
        self.rowconfigure(1, weight=1)

        self.fig, self.ax = plt.subplots(figsize=(9, 3.2))
        self.canvas = FigureCanvasTkAgg(self.fig, master=frame)
        self.canvas.get_tk_widget().pack(fill='both', expand=True, padx=4, pady=4)

    # -----------------------------------------------------------
    # Actions
    # -----------------------------------------------------------

    def reset_defaults(self):
        for key, val in DEFAULTS.items():
            self.vars[key].set(val)
        self.direct_entry_var.set(False)
        self.compute()

    def _read_inputs(self):
        """Parses the entry fields into SI-unit floats. Raises ValueError
        with a field-specific message on bad input, rather than letting
        a bare float() crash propagate to the user as a traceback."""
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
            return

        sigma_o = vals['sigma_o_MPa'] * 1e6
        E = vals['E_GPa'] * 1e9
        R = vals['R_m']
        D = vals['D_m']
        d = vals['d_m']
        t = vals['t_mm'] / 1000.0
        mu = vals['mu']

        direct_entry = self.direct_entry_var.get()

        try:
            result = jt.full_pullin_analysis(R, D, d, t, sigma_o, E, mu,
                                             direct_entry=direct_entry)
        except NotImplementedError as e:
            self.status_var.set(
                'Direct-entry geometry is not implemented in walker1983_jtube.py '
                '(see its docstring for why -- the source paper under-specifies '
                'this case and no reconstruction tried reproduced its stated '
                'result closely enough to trust). Showing the standard '
                'lead-in-section case instead.'
            )
            try:
                result = jt.full_pullin_analysis(R, D, d, t, sigma_o, E, mu,
                                                 direct_entry=False)
            except ValueError as e2:
                self._clear_results()
                self.status_var.set(f'Non-physical inputs: {e2}')
                return
        except ValueError as e:
            self._clear_results()
            self.status_var.set(
                f'Non-physical result for these inputs: {e}\n'
                f'(try a larger R/d or a smaller D/d)'
            )
            return

        self._populate_results(result, vals)
        self._update_plot(result, vals)

    def _clear_results(self):
        for row in self.tree.get_children():
            self.tree.delete(row)

    def _populate_results(self, result, vals):
        self._clear_results()

        s1, s2, s3 = result['stage1'], result['stage2'], result['stage3']
        rows = [
            ('Po (kN)', result['Po'] / 1e3),
            ('Mbar/Mp', s1['Mbar_Mp']),
            ('Stage I: L/d', s1['L_over_d']),
            ('Stage I: l1/d', s1['l1_over_d']),
            ('Stage I: tau (deg)', s1['tau_deg']),
            ('Stage I: P1/Po', s1['P1_over_Po']),
            ('Stage I: P1 (kN)', result['P1'] / 1e3),
            ('Stage II: l2/d', s2['l2_over_d']),
            ('Stage II: P2/Po', s2['P2_over_Po']),
            ('Stage III: Rp/d', s3['Rp_over_d']),
            ("Stage III: P'/Po (peak)", s3['P_prime_over_Po']),
            ('Stage III: P3/Po', s3['P3_over_Po']),
            ('Stage III: P3 (kN)', result['P3'] / 1e3),
        ]

        # is this exactly the Appendix II reference case? (within a tight
        # tolerance on the raw inputs, not the outputs) -- only show a
        # paper comparison column when it's actually meaningful.
        is_reference_case = (
            abs(vals['sigma_o_MPa'] - 400) < 1e-6 and abs(vals['E_GPa'] - 200) < 1e-6
            and abs(vals['R_m'] - 30) < 1e-6 and abs(vals['D_m'] - 0.45) < 1e-6
            and abs(vals['d_m'] - 0.30) < 1e-6 and abs(vals['t_mm'] - 20) < 1e-6
            and abs(vals['mu'] - 0.30) < 1e-6
        )

        for name, value in rows:
            val_str = f'{value:.4g}'
            paper_val = REFERENCE.get(name)
            if is_reference_case and paper_val is not None:
                paper_str = f'{paper_val:.4g}'
                diff_str = f'{100.0*(value-paper_val)/paper_val:+.1f}%'
            else:
                paper_str = '' if not is_reference_case else '-- (not stated)'
                diff_str = ''
            self.tree.insert('', 'end', values=(name, val_str, paper_str, diff_str))

        if not is_reference_case:
            self.status_var.set(
                'Inputs differ from the Appendix II reference case, so no '
                'paper comparison is shown (only the reference case has a '
                'paper-stated value to compare against).'
            )

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

        self.ax.plot(Rds, P1s, color='#1a5c78',
                     label=f'D/d={Dd:.3g}, mu={vals["mu"]:.3g}')

        Rd_current = vals['R_m'] / vals['d_m']
        P1_current = result['stage1']['P1_over_Po'] * 1e3
        self.ax.plot([Rd_current], [P1_current], 'o', color='#97281f',
                     markersize=9, zorder=5,
                     label=f'Current: R/d={Rd_current:.1f}, P1/Po={P1_current:.2f}e-3')

        self.ax.set_xlabel('R/d')
        self.ax.set_ylabel('(P1/Po) x 10^3')
        self.ax.legend(fontsize=8, loc='upper right')
        self.ax.grid(True, alpha=0.3)
        self.fig.tight_layout()
        self.canvas.draw()


def main():
    app = WalkerGUI()
    app.mainloop()


if __name__ == '__main__':
    main()
