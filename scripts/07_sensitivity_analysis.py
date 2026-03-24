#!/usr/bin/env python3
"""
07_sensitivity_analysis.py — P(BH) sensitivity to prior assumptions.

Varies BH mass threshold, M1 prior width, and extinction to show
how the headline P(BH) depends on modelling choices.

Outputs:
  results/sensitivity_results.json
"""

import json, os
import numpy as np
from scipy.optimize import brentq

# ─── Fixed orbital parameters ────────────────────────────────────────
P_DAYS = 44.41
ECC = 0.007
K1 = 80.355  # km/s

N_DRAWS = 200_000
SEED = 42


def mass_function(P, K1, e):
    return 1.0385e-7 * (1 - e**2)**1.5 * K1**3 * P


def run_mc(M1_best, M1_sigma, bh_threshold, fM, seed=SEED):
    """Run MC and return P(BH) for given parameters."""
    rng = np.random.default_rng(seed)
    m1_draws = rng.lognormal(
        mean=np.log(M1_best) - 0.5 * (M1_sigma / M1_best)**2,
        sigma=M1_sigma / M1_best,
        size=N_DRAWS
    )
    cos_i = rng.uniform(0, 1, N_DRAWS)
    sin_i = np.sqrt(1 - cos_i**2)

    m2_draws = []
    for m1, si in zip(m1_draws, sin_i):
        fM_eff = fM / si**3
        try:
            m2 = brentq(lambda m2: m2**3 / (m1 + m2)**2 - fM_eff,
                        0.01, 5000.0)
            m2_draws.append(m2)
        except (ValueError, RuntimeError):
            continue

    m2_arr = np.array(m2_draws)
    p_bh = 100.0 * np.mean(m2_arr > bh_threshold)
    median = float(np.median(m2_arr))
    ci68 = [float(np.percentile(m2_arr, 16)),
            float(np.percentile(m2_arr, 84))]
    return round(p_bh, 1), round(median, 1), [round(c, 1) for c in ci68]


def main():
    fM = mass_function(P_DAYS, K1, ECC)
    print(f'f(M) = {fM:.4f} Msun\n')

    # Fiducial parameters
    M1_FID = 4.40
    SIG_FID = 1.32  # 30%
    BH_FID = 5.0
    AV_FID = 2.85

    # Extinction variations affect M1:
    # Higher AV -> brighter intrinsic -> higher M1
    # Lower AV  -> fainter intrinsic  -> lower M1
    # Approximate: dM1/dAV ~ +0.8 Msun/mag (from HRD slope)
    dm1_per_dav = 0.8

    configs = []

    # 1. Vary BH threshold: 3, 4, 5 Msun
    for bh_thr in [3.0, 4.0, 5.0]:
        label = f'BH > {bh_thr:.0f} Msun, M1=4.40±1.32, AV=2.85'
        p, med, ci = run_mc(M1_FID, SIG_FID, bh_thr, fM)
        configs.append({
            'label': label,
            'BH_threshold': bh_thr,
            'M1': M1_FID,
            'M1_sigma': SIG_FID,
            'A_V': AV_FID,
            'P_BH': p,
            'M2_median': med,
            'M2_68ci': ci,
        })

    # 2. Wider M1 prior (50% uncertainty)
    sig_wide = 2.20  # 50%
    for bh_thr in [3.0, 5.0]:
        label = f'BH > {bh_thr:.0f} Msun, M1=4.40±2.20 (wide), AV=2.85'
        p, med, ci = run_mc(M1_FID, sig_wide, bh_thr, fM)
        configs.append({
            'label': label,
            'BH_threshold': bh_thr,
            'M1': M1_FID,
            'M1_sigma': sig_wide,
            'A_V': AV_FID,
            'P_BH': p,
            'M2_median': med,
            'M2_68ci': ci,
        })

    # 3. Vary extinction: AV = 2.3, 2.85, 3.4
    for av in [2.3, 3.4]:
        m1_shifted = M1_FID + dm1_per_dav * (av - AV_FID)
        sig_shifted = m1_shifted * (SIG_FID / M1_FID)
        label = (f'BH > 5 Msun, M1={m1_shifted:.2f}±{sig_shifted:.2f}, '
                 f'AV={av:.1f}')
        p, med, ci = run_mc(m1_shifted, sig_shifted, 5.0, fM)
        configs.append({
            'label': label,
            'BH_threshold': 5.0,
            'M1': round(m1_shifted, 2),
            'M1_sigma': round(sig_shifted, 2),
            'A_V': av,
            'P_BH': p,
            'M2_median': med,
            'M2_68ci': ci,
        })

    # 4. No B9 anchor: use broader M1 = 3.0 ± 1.5 (generic A/B star)
    label = 'BH > 5 Msun, M1=3.0±1.5 (no B9 anchor), AV=2.85'
    p, med, ci = run_mc(3.0, 1.5, 5.0, fM)
    configs.append({
        'label': label,
        'BH_threshold': 5.0,
        'M1': 3.0,
        'M1_sigma': 1.5,
        'A_V': AV_FID,
        'P_BH': p,
        'M2_median': med,
        'M2_68ci': ci,
    })

    # Print results
    print(f'{"Configuration":<55} {"P(BH)":>7} {"M2 med":>7} {"68% CI":>15}')
    print('-' * 90)
    for c in configs:
        ci_str = f'[{c["M2_68ci"][0]:.1f}, {c["M2_68ci"][1]:.1f}]'
        print(f'{c["label"]:<55} {c["P_BH"]:>6.1f}% '
              f'{c["M2_median"]:>6.1f} {ci_str:>15}')

    # Save
    basedir = os.path.dirname(__file__)
    outpath = os.path.join(basedir, '..', 'results',
                           'sensitivity_results.json')
    with open(outpath, 'w') as f:
        json.dump(configs, f, indent=2)
    print(f'\n✓ Saved: {outpath}')


if __name__ == '__main__':
    main()
