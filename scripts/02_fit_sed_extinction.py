#!/usr/bin/env python3
"""
02_fit_sed_extinction.py — SED analysis with extinction correction for HD 338492.

Uses SIMBAD B9 spectral type to derive interstellar reddening, then
fits a dereddened single-star blackbody to the 10-band photometry
(Gaia G/BP/RP + 2MASS JHKs + WISE W1-W4).

Outputs:
  results/sed_fit_results.json
  paper/figures/fig_sed_extinction.pdf
"""

import json, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

# ─── constants ────────────────────────────────────────────────────────
h_cgs = 6.626e-27     # erg s
c_cgs = 2.998e10      # cm/s
k_cgs = 1.381e-16     # erg/K
Jy = 1e-23            # erg/s/cm²/Hz

# Filter data: zero-point flux (Jy), effective wavelength (μm),
# A(λ)/A_V from Cardelli+1989 / Wang & Chen 2019
FILTERS = {
    'G':  {'fzp': 3228.75, 'lam': 0.622, 'aav': 0.789},
    'BP': {'fzp': 3552.01, 'lam': 0.511, 'aav': 1.002},
    'RP': {'fzp': 2554.95, 'lam': 0.777, 'aav': 0.589},
    'J':  {'fzp': 1594.0,  'lam': 1.235, 'aav': 0.282},
    'H':  {'fzp': 1024.0,  'lam': 1.662, 'aav': 0.175},
    'Ks': {'fzp': 666.7,   'lam': 2.159, 'aav': 0.112},
    'W1': {'fzp': 309.5,   'lam': 3.353, 'aav': 0.065},
    'W2': {'fzp': 171.8,   'lam': 4.603, 'aav': 0.053},
    'W3': {'fzp': 31.67,   'lam': 11.56, 'aav': 0.020},
    'W4': {'fzp': 8.363,   'lam': 22.09, 'aav': 0.010},
}

# ─── HD 338492 data ──────────────────────────────────────────────────
MAGS = {
    'G': 10.04, 'BP': 10.47, 'RP': 9.43,
    'J': 8.384, 'H': 7.813, 'Ks': 7.624,
    'W1': 7.474, 'W2': 7.556, 'W3': 7.487, 'W4': 7.269,
}
BP_RP_OBS = 1.039
SIMBAD_SPT = 'B9'
TEFF_SPT = 10500  # K, from B9 spectral type
TEFF_PHOT = 5376  # K, from observed (reddened) BP-RP colour
DIST_PC = 811.6
DM = 5 * np.log10(DIST_PC / 10)  # distance modulus = 9.55


def blackbody_flux(T, lam_um):
    """Planck function B_ν(T) at wavelength λ in Jy·sr⁻¹-like units."""
    lam_cm = lam_um * 1e-4
    nu = c_cgs / lam_cm
    x = h_cgs * nu / (k_cgs * T)
    if x > 500:
        return 0.0
    return (2 * h_cgs * nu**3 / c_cgs**2) / (np.exp(x) - 1) / Jy


def derive_reddening():
    """Derive E(BP-RP) and A_V from SpT vs photometric Teff mismatch."""
    # Intrinsic (BP-RP)_0 for B9 star ≈ -0.05
    bp_rp_intrinsic = -0.15 + (10500 - TEFF_SPT) / 5000  # ≈ -0.15 for B9
    ebr = max(BP_RP_OBS - bp_rp_intrinsic, 0)
    ag = ebr * 1.89  # A_G / E(BP-RP) ratio from Gaia
    av = ag / 0.789  # A_V = A_G / (A_G/A_V)
    return ebr, ag, av


def fit_sed(av, teff):
    """Fit single-star blackbody SED to dereddened photometry."""
    bands, obs_flux, obs_err, model_flux = [], [], [], []

    for band, mag in MAGS.items():
        f = FILTERS[band]
        # Deredden: subtract extinction in this band
        a_band = f['aav'] * av
        mag_dered = mag - a_band
        flux = f['fzp'] * 10**(-0.4 * mag_dered)
        bands.append(band)
        obs_flux.append(flux)
        obs_err.append(flux * 0.05)  # 5% uncertainty
        model_flux.append(blackbody_flux(teff, f['lam']))

    obs_flux = np.array(obs_flux)
    obs_err = np.array(obs_err)
    model_flux = np.array(model_flux)

    # Least-squares scaling
    scale = np.sum(obs_flux * model_flux / obs_err**2) / \
            np.sum(model_flux**2 / obs_err**2)
    model_scaled = scale * model_flux
    residuals = (obs_flux - model_scaled) / obs_err
    chi2 = np.sum(residuals**2)
    ndof = len(bands) - 2
    chi2_red = chi2 / max(ndof, 1)

    return bands, obs_flux, obs_err, model_scaled, residuals, chi2_red


def make_figure(bands, obs_dered, obs_raw, model_dered, model_raw,
                av, figpath):
    """Create SED figure with and without extinction correction."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    lams = [FILTERS[b]['lam'] for b in bands]

    # Left panel: uncorrected
    ax1.scatter(lams, obs_raw, c='red', s=60, zorder=5, label='Observed')
    ax1.plot(lams, model_raw, 'k--', alpha=0.7, label=f'BB (Teff=10500 K, no extinction)')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel(r'Wavelength ($\mu$m)')
    ax1.set_ylabel('Flux density (Jy)')
    ax1.set_title('Without extinction correction')
    ax1.legend(fontsize=8)
    ax1.xaxis.set_major_formatter(ScalarFormatter())
    ax1.set_xlim(0.3, 30)

    # Right panel: dereddened
    ax2.scatter(lams, obs_dered, c='blue', s=60, zorder=5,
                label=f'Dereddened ($A_V={av:.2f}$)')
    ax2.plot(lams, model_dered, 'k-', alpha=0.7,
             label=f'BB (Teff=10500 K)')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_xlabel(r'Wavelength ($\mu$m)')
    ax2.set_ylabel('Flux density (Jy)')
    ax2.set_title(f'With extinction correction ($A_V={av:.2f}$ mag)')
    ax2.legend(fontsize=8)
    ax2.xaxis.set_major_formatter(ScalarFormatter())
    ax2.set_xlim(0.3, 30)

    fig.suptitle('HD 338492 — SED Analysis', fontsize=13, fontweight='bold')
    fig.tight_layout()
    fig.savefig(figpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  → {figpath}')


def main():
    print('=== SED + Extinction Analysis for HD 338492 ===\n')

    # Derive reddening
    ebr, ag, av = derive_reddening()
    print(f'  Reddening from B9 spectral type:')
    print(f'    E(BP-RP) = {ebr:.3f} mag')
    print(f'    A_G      = {ag:.3f} mag')
    print(f'    A_V      = {av:.3f} mag')

    # Fit dereddened SED
    bands, obs_dered, obs_err, model_dered, res_dered, chi2_dered = \
        fit_sed(av, TEFF_SPT)
    print(f'\n  Dereddened SED fit: χ²_red = {chi2_dered:.2f}')

    # Fit uncorrected SED (for comparison)
    _, obs_raw, _, model_raw, _, chi2_raw = fit_sed(0.0, TEFF_SPT)
    print(f'  Uncorrected SED fit: χ²_red = {chi2_raw:.2f}')
    print(f'  Improvement: {chi2_raw/chi2_dered:.1f}×')

    # Absolute magnitude
    mg_raw = MAGS['G'] - DM
    mg_corr = MAGS['G'] - ag - DM
    print(f'\n  M_G (raw)       = {mg_raw:.2f}')
    print(f'  M_G (corrected) = {mg_corr:.2f}')

    # Save results
    basedir = os.path.dirname(__file__)
    results = {
        'teff_spt': TEFF_SPT,
        'teff_phot': TEFF_PHOT,
        'simbad_sptype': SIMBAD_SPT,
        'E_BP_RP': round(ebr, 4),
        'A_G': round(ag, 4),
        'A_V': round(av, 4),
        'M_G_raw': round(mg_raw, 3),
        'M_G_corrected': round(mg_corr, 3),
        'chi2_red_uncorrected': round(chi2_raw, 2),
        'chi2_red_dereddened': round(chi2_dered, 2),
        'improvement_factor': round(chi2_raw / chi2_dered, 1),
        'bands': bands,
        'residuals_sigma': [round(r, 2) for r in res_dered],
    }
    outpath = os.path.join(basedir, '..', 'results', 'sed_fit_results.json')
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2)
    print(f'\n  ✓ Saved: {outpath}')

    # Make figure
    figpath = os.path.join(basedir, '..', 'paper', 'figures',
                           'fig_sed_extinction.pdf')
    make_figure(bands, obs_dered, obs_raw, model_dered, model_raw,
                av, figpath)

    print('\n=== SED analysis complete ===')


if __name__ == '__main__':
    main()
