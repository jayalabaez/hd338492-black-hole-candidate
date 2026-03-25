#!/usr/bin/env python3
"""
07b_sensitivity_expanded.py — Expanded sensitivity analysis with
spectral-type anchors (B8, B9, A0, A2) and R_V variations.

Adds to the original 07 analysis:
  • B8 anchor (Teff=11400, M1=3.8)
  • A0 anchor (Teff=9700, M1=2.9)
  • A2 anchor (Teff=8840, M1=2.4)
  • R_V = 2.5, 3.5, 4.0 (non-standard extinction)
  • Combined: each anchor × each R_V

Outputs
-------
results/sensitivity_expanded_results.json
"""

import json, os
import numpy as np
from scipy.optimize import brentq

# ─── Fixed orbital parameters ────────────────────────────────────────
P_DAYS = 44.41
ECC    = 0.007
K1     = 80.355   # km/s
N_DRAWS = 50_000
SEED    = 42

# ── Spectral anchors: SpT → (Teff, M1_mean, M1_sigma) ───────────────
# Teff from Pecaut & Mamajek (2013); M1 from Eker+2018 MS relations
# sigma = 30% fractional for catalogue anchors
ANCHORS = {
    'B8':  {'Teff': 11400, 'M1': 3.80, 'sig': 1.14, 'bprp0': -0.18},
    'B9':  {'Teff': 10500, 'M1': 4.40, 'sig': 1.32, 'bprp0': -0.15},
    'A0':  {'Teff':  9700, 'M1': 2.90, 'sig': 0.87, 'bprp0': -0.04},
    'A2':  {'Teff':  8840, 'M1': 2.40, 'sig': 0.72, 'bprp0':  0.09},
}

# ── R_V values & their effect on A_V ────────────────────────────────
# Base: E(BP-RP) = BP_RP_obs - (BP-RP)_0
BP_RP_OBS = 1.039
# A_G / E(BP-RP) ≈ 1.89 for R_V = 3.1 (Jordi+2010)
# For different R_V, the conversion changes.  Approximate:
#   A_G(R_V) ≈ A_G(3.1) × R_V / 3.1  (to first order)
# Then A_V = A_G / 0.789
RV_VALUES = [2.5, 3.1, 3.5, 4.0]

# ── dM1/dAV: how M1 shifts with extinction (from HRD slope) ─────────
DM1_DAV = 0.8  # Msun per mag


def mass_function(P, K1, e):
    """Spectroscopic mass function in solar masses."""
    return 1.0385e-7 * (1 - e**2)**1.5 * K1**3 * P


def compute_av(spt, R_V=3.1):
    """Compute A_V from spectral-type anchor and R_V."""
    bprp0 = ANCHORS[spt]['bprp0']
    ebr = max(BP_RP_OBS - bprp0, 0)
    ag = ebr * 1.89 * (R_V / 3.1)
    av = ag / 0.789
    return round(av, 2)


def run_mc(M1_best, M1_sigma, bh_threshold, fM, seed=SEED):
    """Monte Carlo P(BH) for given parameters."""
    rng = np.random.default_rng(seed)
    mu_ln = np.log(M1_best) - 0.5 * (M1_sigma / M1_best)**2
    sig_ln = M1_sigma / M1_best
    m1_draws = rng.lognormal(mean=mu_ln, sigma=sig_ln, size=N_DRAWS)
    cos_i = rng.uniform(0, 1, N_DRAWS)
    sin_i = np.sqrt(1 - cos_i**2)
    sin_i = np.clip(sin_i, 0.05, 1.0)

    m2_list = []
    for m1, si in zip(m1_draws, sin_i):
        fM_eff = fM / si**3
        if fM_eff > 200:
            continue
        try:
            m2 = brentq(lambda m2: m2**3 / (m1 + m2)**2 - fM_eff,
                        0.01, 500.0, maxiter=200)
            m2_list.append(m2)
        except (ValueError, RuntimeError):
            continue

    m2 = np.array(m2_list)
    p_bh = 100.0 * np.mean(m2 > bh_threshold)
    med = float(np.median(m2))
    ci = [float(np.percentile(m2, 16)), float(np.percentile(m2, 84))]
    # Also compute M2_min at sin i = 1
    m2_min_val = brentq(lambda x: x**3 / (M1_best + x)**2 - fM,
                        0.01, 500.0)
    return {
        'P_BH': round(p_bh, 1),
        'M2_median': round(med, 1),
        'M2_68ci': [round(ci[0], 1), round(ci[1], 1)],
        'M2_min': round(m2_min_val, 2),
        'N_valid': len(m2_list),
    }


def main():
    fM = mass_function(P_DAYS, K1, ECC)
    print(f'f(M) = {fM:.4f} Msun\n')

    all_configs = []

    # ── 1. Each anchor at R_V = 3.1, BH > 5 ─────────────────────────
    print(f'{"Config":<50} {"AV":>5} {"M2min":>6} {"P(BH)":>7} '
          f'{"Med":>5} {"68% CI":>15}')
    print('-' * 95)

    for spt, info in ANCHORS.items():
        av = compute_av(spt, 3.1)
        m1 = info['M1']
        sig = info['sig']
        r = run_mc(m1, sig, 5.0, fM)
        label = f'{spt} anchor, M1={m1}±{sig}, R_V=3.1, AV={av}'
        ci_str = f'[{r["M2_68ci"][0]}, {r["M2_68ci"][1]}]'
        print(f'  {label:<48} {av:>5.2f} {r["M2_min"]:>6.2f} '
              f'{r["P_BH"]:>6.1f}% {r["M2_median"]:>5.1f} {ci_str:>15}')
        all_configs.append({
            'label': label, 'SpT': spt, 'Teff': info['Teff'],
            'M1': m1, 'M1_sigma': sig, 'R_V': 3.1, 'A_V': av,
            'BH_threshold': 5.0, **r
        })

    # ── 2. B9 anchor with BH thresholds 3, 4, 5 ─────────────────────
    print()
    for bh_thr in [3.0, 4.0]:
        av = compute_av('B9', 3.1)
        m1, sig = 4.40, 1.32
        r = run_mc(m1, sig, bh_thr, fM)
        label = f'B9, BH>{bh_thr}, M1={m1}±{sig}'
        ci_str = f'[{r["M2_68ci"][0]}, {r["M2_68ci"][1]}]'
        print(f'  {label:<48} {av:>5.2f} {r["M2_min"]:>6.2f} '
              f'{r["P_BH"]:>6.1f}% {r["M2_median"]:>5.1f} {ci_str:>15}')
        all_configs.append({
            'label': label, 'SpT': 'B9', 'Teff': 10500,
            'M1': m1, 'M1_sigma': sig, 'R_V': 3.1, 'A_V': av,
            'BH_threshold': bh_thr, **r
        })

    # ── 3. R_V variations (B9 anchor) ────────────────────────────────
    print()
    for rv in RV_VALUES:
        if abs(rv - 3.1) < 0.01:
            continue  # already done
        av = compute_av('B9', rv)
        # Adjust M1 for changed AV
        m1_adj = 4.40 + DM1_DAV * (av - compute_av('B9', 3.1))
        sig_adj = m1_adj * 0.30
        r = run_mc(m1_adj, sig_adj, 5.0, fM)
        label = f'B9, R_V={rv}, AV={av}, M1={m1_adj:.2f}±{sig_adj:.2f}'
        ci_str = f'[{r["M2_68ci"][0]}, {r["M2_68ci"][1]}]'
        print(f'  {label:<48} {av:>5.2f} {r["M2_min"]:>6.2f} '
              f'{r["P_BH"]:>6.1f}% {r["M2_median"]:>5.1f} {ci_str:>15}')
        all_configs.append({
            'label': label, 'SpT': 'B9', 'Teff': 10500,
            'M1': round(m1_adj, 2), 'M1_sigma': round(sig_adj, 2),
            'R_V': rv, 'A_V': av, 'BH_threshold': 5.0, **r
        })

    # ── 4. Cross: each anchor × each non-standard R_V ───────────────
    print('\n  Cross-validation (anchor × R_V):')
    for spt in ['B8', 'A0']:
        for rv in [2.5, 4.0]:
            av = compute_av(spt, rv)
            base_m1 = ANCHORS[spt]['M1']
            m1_adj = base_m1 + DM1_DAV * (av - compute_av(spt, 3.1))
            sig_adj = m1_adj * 0.30
            r = run_mc(max(m1_adj, 0.5), sig_adj, 5.0, fM)
            label = f'{spt}, R_V={rv}, AV={av}, M1={m1_adj:.2f}'
            ci_str = f'[{r["M2_68ci"][0]}, {r["M2_68ci"][1]}]'
            print(f'  {label:<48} {av:>5.2f} {r["M2_min"]:>6.2f} '
                  f'{r["P_BH"]:>6.1f}% {r["M2_median"]:>5.1f} {ci_str:>15}')
            all_configs.append({
                'label': label, 'SpT': spt, 'Teff': ANCHORS[spt]['Teff'],
                'M1': round(m1_adj, 2), 'M1_sigma': round(sig_adj, 2),
                'R_V': rv, 'A_V': av, 'BH_threshold': 5.0, **r
            })

    # ── 5. Extreme: no SpT anchor at all ─────────────────────────────
    r = run_mc(3.0, 1.5, 5.0, fM)
    label = 'No anchor, M1=3.0±1.5'
    all_configs.append({
        'label': label, 'SpT': 'none', 'Teff': None,
        'M1': 3.0, 'M1_sigma': 1.5, 'R_V': 3.1, 'A_V': 2.85,
        'BH_threshold': 5.0, **r
    })
    print(f'\n  {label}: P(BH)={r["P_BH"]}%, M2_min={r["M2_min"]}')

    # ── Summary ──────────────────────────────────────────────────────
    pbh_vals = [c['P_BH'] for c in all_configs]
    print(f'\n  P(BH) range: {min(pbh_vals):.1f}% – {max(pbh_vals):.1f}%')
    print(f'  Total configurations tested: {len(all_configs)}')

    # ── Save ─────────────────────────────────────────────────────────
    basedir = os.path.dirname(__file__)
    outpath = os.path.join(basedir, '..', 'results',
                           'sensitivity_expanded_results.json')
    with open(outpath, 'w') as f:
        json.dump(all_configs, f, indent=2)
    print(f'\n  ✓ Saved: {outpath}')


if __name__ == '__main__':
    main()
