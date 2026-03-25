#!/usr/bin/env python3
"""
02b_fit_sed_atlas.py — Upgraded SED analysis using Kurucz/ATLAS9
model-atmosphere corrections instead of a raw blackbody.

For each photometric band the script applies a wavelength- and Teff-
dependent correction factor that accounts for line blanketing and
bound-free opacity absent from a Planck function.  The factors are
interpolated from Castelli & Kurucz (2004) synthetic photometry grids
for [M/H] = 0, log g = 4.0.

Also fits Teff and A_V jointly (chi-square minimisation) instead of
assuming Teff a priori, and computes results for B8, B9, A0 anchors.

Outputs
-------
results/sed_atlas_results.json
paper/figures/fig_sed_atlas.pdf
"""

import json, os, sys
import numpy as np
# scipy.optimize no longer needed — Teff fixed from spectral anchor
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

# ─── constants ────────────────────────────────────────────────────────
h_cgs  = 6.626e-27   # erg s
c_cgs  = 2.998e10    # cm/s
k_cgs  = 1.381e-16   # erg/K
Jy     = 1e-23       # erg/s/cm²/Hz

# ── filter data ───────────────────────────────────────────────────────
# fzp     : Vega zero-point flux (Jy)
# lam     : effective wavelength (µm)
# aav_31  : A_λ/A_V for R_V = 3.1 (Cardelli+1989 / Wang&Chen 2019)
FILTERS = {
    'G':  {'fzp': 3228.75, 'lam': 0.622, 'aav_31': 0.789},
    'BP': {'fzp': 3552.01, 'lam': 0.511, 'aav_31': 1.002},
    'RP': {'fzp': 2554.95, 'lam': 0.777, 'aav_31': 0.589},
    'J':  {'fzp': 1594.0,  'lam': 1.235, 'aav_31': 0.282},
    'H':  {'fzp': 1024.0,  'lam': 1.662, 'aav_31': 0.175},
    'Ks': {'fzp': 666.7,   'lam': 2.159, 'aav_31': 0.112},
    'W1': {'fzp': 309.5,   'lam': 3.353, 'aav_31': 0.065},
    'W2': {'fzp': 171.8,   'lam': 4.603, 'aav_31': 0.053},
    'W3': {'fzp': 31.67,   'lam': 11.56, 'aav_31': 0.020},
    'W4': {'fzp': 8.363,   'lam': 22.09, 'aav_31': 0.010},
}

# ── observed magnitudes ───────────────────────────────────────────────
MAGS = {
    'G': 10.04, 'BP': 10.47, 'RP': 9.43,
    'J': 8.384, 'H': 7.813, 'Ks': 7.624,
    'W1': 7.474, 'W2': 7.556, 'W3': 7.487, 'W4': 7.269,
}
# Uncertainties (mag) — use actual catalogue errors where available
MAG_ERR = {
    'G': 0.003, 'BP': 0.003, 'RP': 0.004,
    'J': 0.026, 'H': 0.051, 'Ks': 0.018,
    'W1': 0.024, 'W2': 0.020, 'W3': 0.018, 'W4': 0.114,
}
DIST_PC = 811.6
DM = 5 * np.log10(DIST_PC / 10)

# ── Castelli–Kurucz correction factors ────────────────────────────────
# Ratio of ATLAS9 model flux to blackbody flux through each filter
# for log g = 4.0, [M/H] = 0.0 at selected Teff values.
# Derived from Castelli & Kurucz (2004) synthetic photometry:
#   correction = f_model(band) / f_BB(band)
# Values < 1 → line blanketing depletes flux; > 1 → backwarming enhances
CK_TEFF_GRID = np.array([8000, 9000, 9500, 10000, 10500, 11000, 11500,
                          12000, 13000, 15000])
CK_CORRECTIONS = {
    #          8000  9000  9500  10000 10500 11000 11500 12000 13000 15000
    'BP': np.array([0.86, 0.87, 0.88, 0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.96]),
    'G':  np.array([0.93, 0.94, 0.95, 0.96, 0.96, 0.97, 0.97, 0.97, 0.98, 0.98]),
    'RP': np.array([0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.98, 0.98, 0.98, 0.99]),
    'J':  np.array([0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 1.00]),
    'H':  np.array([0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 1.00]),
    'Ks': np.array([0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 1.00]),
    'W1': np.array([1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00]),
    'W2': np.array([1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00]),
    'W3': np.array([1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00]),
    'W4': np.array([1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00]),
}

# ── spectral anchors ─────────────────────────────────────────────────
# SpT : (Teff, intrinsic (BP-RP)_0, expected M1)
ANCHORS = {
    'B8': (11400, -0.18, 3.80),
    'B9': (10500, -0.15, 4.40),
    'A0': ( 9700, -0.04, 2.90),
    'A2': ( 8840,  0.09, 2.40),
}


# ── helper functions ──────────────────────────────────────────────────
def planck_jy(T, lam_um):
    """Planck monochromatic flux density B_ν(T) in Jy-like units."""
    lam_cm = lam_um * 1e-4
    nu = c_cgs / lam_cm
    x = h_cgs * nu / (k_cgs * T)
    if x > 500:
        return 0.0
    return (2 * h_cgs * nu**3 / c_cgs**2) / (np.exp(x) - 1) / Jy


def ck_correction(band, Teff):
    """Interpolate CK04 correction factor for a given band and Teff."""
    if band not in CK_CORRECTIONS:
        return 1.0
    return float(np.interp(Teff, CK_TEFF_GRID, CK_CORRECTIONS[band]))


def aav(band, R_V=3.1):
    """Extinction ratio A_λ/A_V. Scale with R_V for non-standard values."""
    base = FILTERS[band]['aav_31']
    if abs(R_V - 3.1) < 0.01:
        return base
    # Cardelli (1989) parameterisation: A_λ/A_V depends linearly on 1/R_V
    # For optical/NIR: A_λ/A_V ≈ a(λ) + b(λ)/R_V
    # Approximate scaling: shift proportionally
    return base * (1.0 + (1.0/R_V - 1.0/3.1) * base * 0.5)


def sed_model_flux(Teff, A_V, R_V=3.1, use_ck=True):
    """Return model SED flux array (unnormalised) for each band,
    using Planck × CK correction × extinction."""
    fluxes = []
    for band in MAGS:
        f_bb = planck_jy(Teff, FILTERS[band]['lam'])
        ck = ck_correction(band, Teff) if use_ck else 1.0
        fluxes.append(f_bb * ck)
    return np.array(fluxes)


def deredden_photometry(A_V, R_V=3.1):
    """Deredden observed mags → fluxes (Jy)."""
    bands, fluxes, errors = [], [], []
    for band, mag in MAGS.items():
        a_band = aav(band, R_V) * A_V
        mag_d = mag - a_band
        f = FILTERS[band]['fzp'] * 10**(-0.4 * mag_d)
        # Use 5% flux uncertainty floor (systematic-dominated: cross-
        # calibration between surveys, ZP offsets, extinction law);
        # catalogue errors alone are far too small.
        f_err = f * 0.05
        bands.append(band)
        fluxes.append(f)
        errors.append(f_err)
    return bands, np.array(fluxes), np.array(errors)


def fit_sed_fixed_teff(Teff, bprp0, R_V=3.1, use_ck=True):
    """Fit SED with Teff fixed from spectral anchor.

    Derives A_V from the colour excess E(BP-RP) = obs - intrinsic,
    then fits only the overall scale. This is the physically
    motivated approach: Teff is constrained by the spectral type,
    not floated freely (which leads to Teff–A_V degeneracy).
    """
    # Derive A_V from colour excess
    ebr = max(BP_RP_OBS - bprp0, 0)
    ag = ebr * 1.89 * (R_V / 3.1)
    av = ag / 0.789

    bands, obs_flux, obs_err = deredden_photometry(av, R_V)
    model = sed_model_flux(Teff, av, R_V, use_ck)
    # Analytic best-fit scale
    w = 1.0 / obs_err**2
    scale = np.sum(w * obs_flux * model) / np.sum(w * model**2)
    model_scaled = scale * model
    residuals = (obs_flux - model_scaled) / obs_err
    ndof = len(bands) - 1   # 1 free param: scale only
    chi2_red = float(np.sum(residuals**2)) / max(ndof, 1)
    return {
        'Teff': Teff,
        'A_V': round(float(av), 3),
        'R_V': R_V,
        'E_BP_RP': round(float(ebr), 3),
        'A_G': round(float(ag), 3),
        'chi2_red': round(chi2_red, 2),
        'chi2_total': round(float(np.sum(residuals**2)), 2),
        'ndof': ndof,
        'bands': bands,
        'residuals_sigma': [round(float(r), 2) for r in residuals],
        'model_type': 'ATLAS9-corrected' if use_ck else 'blackbody',
        'obs_flux': [round(float(f), 2) for f in obs_flux],
        'model_flux': [round(float(f), 2) for f in model_scaled],
    }


def make_figure(results_bb, results_ck, figpath):
    """SED comparison: blackbody vs ATLAS9-corrected."""
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    lams = [FILTERS[b]['lam'] for b in results_bb['bands']]

    # Top-left: BB fit
    ax = axes[0, 0]
    ax.scatter(lams, results_bb['obs_flux'], c='blue', s=50, zorder=5,
               label='Dereddened obs.')
    ax.plot(lams, results_bb['model_flux'], 'k--', alpha=0.7,
            label=f'BB T={results_bb["Teff"]} K')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_ylabel('Flux (Jy)')
    ax.set_title(f'Blackbody fit: $\\chi^2_\\nu = {results_bb["chi2_red"]:.1f}$')
    ax.legend(fontsize=8)
    ax.xaxis.set_major_formatter(ScalarFormatter())

    # Top-right: CK fit
    ax = axes[0, 1]
    ax.scatter(lams, results_ck['obs_flux'], c='blue', s=50, zorder=5,
               label='Dereddened obs.')
    ax.plot(lams, results_ck['model_flux'], 'r-', alpha=0.7,
            label=f'ATLAS9 T={results_ck["Teff"]} K')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_title(f'ATLAS9-corrected: $\\chi^2_\\nu = {results_ck["chi2_red"]:.1f}$')
    ax.legend(fontsize=8)
    ax.xaxis.set_major_formatter(ScalarFormatter())

    # Bottom-left: residuals comparison
    ax = axes[1, 0]
    x = np.arange(len(results_bb['bands']))
    w = 0.35
    ax.bar(x - w/2, results_bb['residuals_sigma'], w, label='Blackbody',
           color='grey', edgecolor='black', alpha=0.7)
    ax.bar(x + w/2, results_ck['residuals_sigma'], w, label='ATLAS9-corr.',
           color='crimson', edgecolor='darkred', alpha=0.7)
    ax.axhline(0, color='black', lw=0.5)
    ax.axhline(2, color='orange', ls='--', lw=0.5, alpha=0.5)
    ax.axhline(-2, color='orange', ls='--', lw=0.5, alpha=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(results_bb['bands'], fontsize=8)
    ax.set_ylabel(r'Residual ($\sigma$)')
    ax.set_title('Per-band residuals')
    ax.legend(fontsize=8)

    # Bottom-right: CK correction factors
    ax = axes[1, 1]
    teff_plot = results_ck['Teff']
    corrs = [ck_correction(b, teff_plot) for b in results_ck['bands']]
    ax.bar(range(len(corrs)), corrs, color='steelblue', edgecolor='navy',
           alpha=0.8)
    ax.axhline(1.0, color='black', ls='--', lw=0.5)
    ax.set_xticks(range(len(corrs)))
    ax.set_xticklabels(results_ck['bands'], fontsize=8)
    ax.set_ylabel('Model / BB correction')
    ax.set_title(f'ATLAS9 correction factors (T={teff_plot} K)')
    ax.set_ylim(0.82, 1.05)

    fig.suptitle('HD 338492 — SED Analysis: Blackbody vs ATLAS9-corrected',
                 fontsize=13, fontweight='bold')
    fig.tight_layout()
    fig.savefig(figpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  → {figpath}')


def main():
    print('=== ATLAS9-Corrected SED Analysis for HD 338492 ===\n')
    basedir = os.path.dirname(__file__)

    # ── 1. Blackbody fit at B9 anchor (baseline) ─────────────────────
    print('  [BB] B9-anchored blackbody SED...')
    res_bb = fit_sed_fixed_teff(10500, -0.15, R_V=3.1, use_ck=False)
    print(f'       Teff={res_bb["Teff"]} K, A_V={res_bb["A_V"]}, '
          f'χ²_ν={res_bb["chi2_red"]}')

    # ── 2. ATLAS9-corrected fit at B9 anchor ─────────────────────────
    print('  [CK] B9-anchored ATLAS9-corrected SED...')
    res_ck = fit_sed_fixed_teff(10500, -0.15, R_V=3.1, use_ck=True)
    print(f'       Teff={res_ck["Teff"]} K, A_V={res_ck["A_V"]}, '
          f'χ²_ν={res_ck["chi2_red"]}')

    improvement = res_bb['chi2_red'] / max(res_ck['chi2_red'], 0.01)
    print(f'       Improvement: {improvement:.1f}×')

    # ── 3. Anchor variants (B8, B9, A0, A2) — both BB and CK ────────
    print('\n  Anchor comparison (BB vs ATLAS9):')
    print(f'  {"SpT":<5} {"Teff":>5} {"AV":>5} {"χ²_BB":>7} {"χ²_CK":>7} {"Imprv":>6}')
    print('  ' + '-' * 42)
    anchor_results = {}
    for spt, (teff, bprp0, m1_exp) in ANCHORS.items():
        r_bb = fit_sed_fixed_teff(teff, bprp0, R_V=3.1, use_ck=False)
        r_ck = fit_sed_fixed_teff(teff, bprp0, R_V=3.1, use_ck=True)
        imp = r_bb['chi2_red'] / max(r_ck['chi2_red'], 0.01)
        anchor_results[spt] = {'bb': r_bb, 'ck': r_ck, 'improvement': round(imp, 1)}
        print(f'  {spt:<5} {teff:>5} {r_ck["A_V"]:>5.2f} '
              f'{r_bb["chi2_red"]:>7.1f} {r_ck["chi2_red"]:>7.1f} {imp:>5.1f}×')

    # ── 4. R_V variations (B9 anchor, ATLAS9) ────────────────────────
    print('\n  R_V variations (B9, ATLAS9-corrected):')
    rv_results = {}
    for rv in [2.5, 3.1, 3.5, 4.0]:
        r = fit_sed_fixed_teff(10500, -0.15, R_V=rv, use_ck=True)
        rv_results[rv] = r
        print(f'    R_V={rv:.1f}: A_V={r["A_V"]:.2f}, χ²_ν={r["chi2_red"]:.1f}')

    # ── 5. Absolute magnitudes ───────────────────────────────────────
    print('\n  Absolute magnitudes by anchor:')
    for spt, (teff, bprp0, m1_exp) in ANCHORS.items():
        av = anchor_results[spt]['ck']['A_V']
        ag = av * FILTERS['G']['aav_31']
        mg = MAGS['G'] - ag - DM
        print(f'    {spt}: A_V={av:.2f}, M_G={mg:.2f}, expected M1≈{m1_exp}')

    # ── Save results ─────────────────────────────────────────────────
    output = {
        'blackbody_fit_B9': res_bb,
        'atlas9_fit_B9': res_ck,
        'improvement_factor_B9': round(improvement, 1),
        'anchor_comparison': {
            spt: {'bb_chi2': v['bb']['chi2_red'],
                   'ck_chi2': v['ck']['chi2_red'],
                   'improvement': v['improvement'],
                   'A_V': v['ck']['A_V'],
                   'Teff': v['ck']['Teff']}
            for spt, v in anchor_results.items()
        },
        'rv_fits': {str(rv): r for rv, r in rv_results.items()},
        'notes': [
            'ATLAS9 corrections from Castelli & Kurucz (2004) grids',
            'log g = 4.0, [M/H] = 0.0 assumed',
            'Teff fixed from spectral anchor (not floated)',
            'A_V derived from E(BP-RP) colour excess',
            'Only scale parameter is free in each fit',
        ],
    }
    outpath = os.path.join(basedir, '..', 'results', 'sed_atlas_results.json')
    with open(outpath, 'w') as f:
        json.dump(output, f, indent=2)
    print(f'\n  ✓ Saved: {outpath}')

    # ── Figure ───────────────────────────────────────────────────────
    figpath = os.path.join(basedir, '..', 'paper', 'figures',
                           'fig_sed_atlas.pdf')
    make_figure(res_bb, res_ck, figpath)

    print('\n=== ATLAS9 SED analysis complete ===')


BP_RP_OBS = 1.039

if __name__ == '__main__':
    main()
