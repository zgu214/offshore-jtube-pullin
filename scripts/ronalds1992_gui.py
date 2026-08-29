# -*- coding: utf-8 -*-
"""
Desktop GUI for ronalds1992_jtube.py.

Same convention as walker1983_gui.py: tkinter (stdlib) + an embedded
matplotlib plot, no extra install beyond matplotlib/numpy. Two tabs,
matching the two independent calculation modes ronalds1992_jtube.py
provides:

  Tab 1 -- Multi-bend tension growth (eq. 1-3, 10-11): back tension,
           friction coefficient, and a list of bend angles in;
           tension after each bend out, plus a plot of tension vs.
           cumulative angle.
  Tab 2 -- Single-bend support response (eq. 5-9, Fig. 5): the
           nondimensional lengths gamma/lambda/xi and a support-type
           selector (Free / Guide / Full support) in; eta, MB, the
           angle and value of the max span moment out.

Unlike walker1983_gui.py there is no paper-stated worked example to
compare against (see ronalds1992_jtube.py's module docstring) --
Tab 2 instead always shows the eq. (6) zero-length reference value
(eta=2/pi) alongside the current GUIDE-case result, since that is the
one number this module's guide-support formula is actually checked
against.

Run:
    python ronalds1992_gui.py
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
import ronalds1992_jtube as rn


DEFAULT_BENDS = [90.0, 8.5, 8.5]   # matches the Fan et al. (2013) multi-bend example


class RonaldsGUI(tk.Tk):

    def __init__(self):
        super().__init__()
        self.title('Ronalds (1992) -- J-Tube Design for Flexible Umbilicals')
        self.geometry('1080x680')
        self.minsize(920, 600)
        self.protocol("WM_DELETE_WINDOW", self._on_close)

        notebook = ttk.Notebook(self)
        notebook.pack(fill='both', expand=True, padx=6, pady=6)

        self.tab_multibend = ttk.Frame(notebook)
        self.tab_support = ttk.Frame(notebook)
        notebook.add(self.tab_multibend, text='Multi-bend tension growth')
        notebook.add(self.tab_support, text='Single-bend support response')

        self._build_multibend_tab()
        self._build_support_tab()

        self.compute_multibend()
        self.compute_support()

    def _on_close(self):
        """
        Explicitly close both embedded matplotlib figures and force the
        process to exit. Without this, closing the window can leave the
        Python process running in the background -- pyplot keeps its own
        global reference to each figure created via plt.subplots(), and
        plain self.destroy() is not reliably enough to let the
        interpreter shut down on its own.
        """
        plt.close('all')
        self.destroy()
        os._exit(0)

    # =================================================================
    # Tab 1: multi-bend tension growth
    # =================================================================

    def _build_multibend_tab(self):
        tab = self.tab_multibend
        tab.columnconfigure(1, weight=1)
        tab.rowconfigure(2, weight=1)

        inputs = ttk.LabelFrame(tab, text='Inputs')
        inputs.grid(row=0, column=0, sticky='nsew', padx=8, pady=8)

        ttk.Label(inputs, text='Back tension T0 (kN)').grid(row=0, column=0, sticky='w', padx=4, pady=3)
        self.T0_var = tk.StringVar(value='100')
        ttk.Entry(inputs, textvariable=self.T0_var, width=10).grid(row=0, column=1, padx=4, pady=3)

        ttk.Label(inputs, text='Friction coefficient, mu').grid(row=1, column=0, sticky='w', padx=4, pady=3)
        self.mu_mb_var = tk.StringVar(value='0.30')
        ttk.Entry(inputs, textvariable=self.mu_mb_var, width=10).grid(row=1, column=1, padx=4, pady=3)

        ttk.Label(inputs, text='Bend angles (deg), one per line,\n'
                              'in pull-in order (nearest back\n'
                              'tension first):').grid(
            row=2, column=0, columnspan=2, sticky='w', padx=4, pady=(10, 2))
        self.bends_text = tk.Text(inputs, width=14, height=8)
        self.bends_text.grid(row=3, column=0, columnspan=2, padx=4, pady=2)
        self.bends_text.insert('1.0', '\n'.join(str(a) for a in DEFAULT_BENDS))

        ttk.Button(inputs, text='Compute', command=self.compute_multibend).grid(
            row=4, column=0, columnspan=2, pady=(10, 4), sticky='ew')
        ttk.Button(inputs, text='Reset to Fan et al. (2013) example geometry',
                  command=self._reset_multibend_defaults).grid(
            row=5, column=0, columnspan=2, pady=2, sticky='ew')

        self.mb_status_var = tk.StringVar(value='')
        ttk.Label(inputs, textvariable=self.mb_status_var, foreground='#97281f',
                 wraplength=180, justify='left').grid(
            row=6, column=0, columnspan=2, sticky='w', padx=4, pady=(8, 4))

        results = ttk.LabelFrame(tab, text='Tension after each bend')
        results.grid(row=0, column=1, sticky='nsew', padx=8, pady=8)
        columns = ('bend', 'angle_deg', 'tension_kN', 'ratio')
        self.mb_tree = ttk.Treeview(results, columns=columns, show='headings', height=10)
        headers = {'bend': 'Bend #', 'angle_deg': 'Angle (deg)',
                  'tension_kN': 'Tension (kN)', 'ratio': 'T / T0'}
        for c in columns:
            self.mb_tree.heading(c, text=headers[c])
            self.mb_tree.column(c, width=110, anchor='center')
        self.mb_tree.pack(fill='both', expand=True, padx=4, pady=4)

        plot_frame = ttk.LabelFrame(tab, text='Tension vs. cumulative bend angle')
        plot_frame.grid(row=1, column=0, columnspan=2, sticky='nsew', padx=8, pady=(0, 8))
        tab.rowconfigure(1, weight=1)
        self.mb_fig, self.mb_ax = plt.subplots(figsize=(8, 2.6))
        self.mb_canvas = FigureCanvasTkAgg(self.mb_fig, master=plot_frame)
        self.mb_canvas.get_tk_widget().pack(fill='both', expand=True, padx=4, pady=4)

    def _reset_multibend_defaults(self):
        self.T0_var.set('100')
        self.mu_mb_var.set('0.30')
        self.bends_text.delete('1.0', 'end')
        self.bends_text.insert('1.0', '\n'.join(str(a) for a in DEFAULT_BENDS))
        self.compute_multibend()

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
            self.mb_tree.insert('', 'end', values=(
                i, f'{a_deg:.2f}', f'{T/1e3:.3f}', f'{T/T0:.4f}'))

        self.mb_ax.clear()
        cum_angle = np.concatenate([[0.0], np.cumsum(angles_deg)])
        self.mb_ax.plot(cum_angle, np.array(tensions) / 1e3, marker='o', color='#1a5c78')
        self.mb_ax.set_xlabel('Cumulative bend angle (deg)')
        self.mb_ax.set_ylabel('Tension (kN)')
        self.mb_ax.grid(True, alpha=0.3)
        self.mb_fig.tight_layout()
        self.mb_canvas.draw()

    # =================================================================
    # Tab 2: single-bend support response
    # =================================================================

    def _build_support_tab(self):
        tab = self.tab_support
        tab.columnconfigure(1, weight=1)
        tab.rowconfigure(1, weight=1)

        inputs = ttk.LabelFrame(tab, text='Inputs')
        inputs.grid(row=0, column=0, sticky='nsew', padx=8, pady=8)

        ttk.Label(inputs, text='Support type at A').grid(row=0, column=0, sticky='w', padx=4, pady=3)
        self.support_var = tk.StringVar(value='Guide')
        ttk.Combobox(inputs, textvariable=self.support_var, state='readonly',
                    values=['Free', 'Guide', 'Full support'], width=12).grid(
            row=0, column=1, padx=4, pady=3)
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
        ttk.Button(inputs, text='Reset to zero-length (eq. 6) case',
                  command=self._reset_support_defaults).grid(
            row=5, column=0, columnspan=2, pady=2, sticky='ew')

        self.sup_status_var = tk.StringVar(value='')
        ttk.Label(inputs, textvariable=self.sup_status_var, foreground='#97281f',
                 wraplength=180, justify='left').grid(
            row=6, column=0, columnspan=2, sticky='w', padx=4, pady=(8, 4))

        note = ('Note: gamma/lambda/xi only affect the GUIDE case here.\n'
               'FREE and FULL SUPPORT results below do not depend on them\n'
               '(FULL SUPPORT is unconditionally V=0, H=TL, MB=0 -- see\n'
               'ronalds1992_jtube.py). Deflection formulas are not\n'
               'implemented -- see the module docstring for why.')
        ttk.Label(inputs, text=note, foreground='#697680',
                 wraplength=200, justify='left', font=('TkDefaultFont', 8)).grid(
            row=7, column=0, columnspan=2, sticky='w', padx=4, pady=(10, 4))

        results = ttk.LabelFrame(tab, text='Result  (normalised by TL and r, where applicable)')
        results.grid(row=0, column=1, sticky='nsew', padx=8, pady=8)
        columns = ('quantity', 'value', 'reference')
        self.sup_tree = ttk.Treeview(results, columns=columns, show='headings', height=10)
        headers = {'quantity': 'Quantity', 'value': 'Value',
                  'reference': "Paper's eq. (6) reference"}
        widths = {'quantity': 200, 'value': 110, 'reference': 200}
        for c in columns:
            self.sup_tree.heading(c, text=headers[c])
            self.sup_tree.column(c, width=widths[c], anchor='w' if c != 'value' else 'center')
        self.sup_tree.pack(fill='both', expand=True, padx=4, pady=4)

        plot_frame = ttk.LabelFrame(tab, text='eta = V/TL vs. lambda  (at the current gamma, xi)')
        plot_frame.grid(row=1, column=0, columnspan=2, sticky='nsew', padx=8, pady=(0, 8))
        self.sup_fig, self.sup_ax = plt.subplots(figsize=(8, 2.6))
        self.sup_canvas = FigureCanvasTkAgg(self.sup_fig, master=plot_frame)
        self.sup_canvas.get_tk_widget().pack(fill='both', expand=True, padx=4, pady=4)

    def _reset_support_defaults(self):
        self.gamma_var.set('0.0')
        self.lam_var.set('0.0')
        self.xi_var.set('0.0')
        self.support_var.set('Guide')
        self.compute_support()

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
            self.sup_status_var.set('gamma, lambda, xi must all be >= 0 '
                                    "(the paper's own stated restriction).")
            return

        support = self.support_var.get()

        for row in self.sup_tree.get_children():
            self.sup_tree.delete(row)

        if support == 'Full support':
            r = rn.full_support_response()
            rows = [
                ('eta = V/TL', r['eta'], ''),
                ('MB/(TL*r)', r['MB_over_TLr'], ''),
                ('H/TL', r['H_over_TL'], ''),
            ]
        elif support == 'Free':
            MB = rn.free_support_MB_over_TLr(xi)
            rows = [
                ('MB/(TL*r)', MB, ''),
                ('eta = V/TL', 0.0, '(V=0 for the free case)'),
                ('H/TL', 0.0, '(H=0 for the free case)'),
            ]
        else:  # Guide
            r = rn.guide_support_response(gamma, lam, xi)
            ref_note = ''
            if gamma == 0.0 and lam == 0.0 and xi == 0.0:
                ref_note = f"matches eq.(6): 2/pi = {2/math.pi:.6f}"
            rows = [
                ('eta = V/TL', r['eta'], ref_note),
                ('MB/(TL*r)', r['MB_over_TLr'], ''),
                ('theta at max Ms (deg)', r['theta_max_deg'], ''),
                ('MS/(TL*r)', r['MS_over_TLr'], ''),
            ]

        for name, value, note in rows:
            self.sup_tree.insert('', 'end', values=(name, f'{value:.5f}', note))

        self._update_support_plot(gamma, xi, support)

    def _update_support_plot(self, gamma, xi, support):
        self.sup_ax.clear()

        if support == 'Guide':
            lam_range = np.linspace(0.0, 8.0, 81)
            etas = [rn.guide_support_eta(gamma, lam, xi) for lam in lam_range]
            self.sup_ax.plot(lam_range, etas, color='#1a5c78', label='eta(lambda)')
            lam_current = max(0.0, float(self.lam_var.get().strip() or 0.0))
            eta_current = rn.guide_support_eta(gamma, lam_current, xi)
            self.sup_ax.plot([lam_current], [eta_current], 'o', color='#97281f',
                             markersize=8, zorder=5, label='current')
            self.sup_ax.axhline(2 / math.pi, color='#8f5307', linestyle=':',
                                alpha=0.7, label='eq.(6): 2/pi (gamma=lambda=xi=0)')
            self.sup_ax.set_ylabel('eta = V/TL')
        else:
            self.sup_ax.text(0.5, 0.5, f'"{support}" support: eta is fixed '
                                        f'(not a function of lambda)\nsee the '
                                        f'results table',
                             ha='center', va='center', transform=self.sup_ax.transAxes,
                             fontsize=10, color='#697680')

        self.sup_ax.set_xlabel('lambda')
        self.sup_ax.grid(True, alpha=0.3)
        self.sup_ax.legend(fontsize=8) if support == 'Guide' else None
        self.sup_fig.tight_layout()
        self.sup_canvas.draw()


def main():
    app = RonaldsGUI()
    app.mainloop()


if __name__ == '__main__':
    main()
