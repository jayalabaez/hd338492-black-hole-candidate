#!/usr/bin/env python3
"""
03_compute_mass_posterior.py — Bayesian mass posterior for HD 338492.

Computes minimum companion mass from the spectroscopic mass function,
then performs Monte Carlo inclination debiasing with M1 uncertainty
propagation (2×10⁵ draws).

Outputs:
  results/mass_posterior_results.json
  paper/figures/fig_mass_posterior.pdf
"""

import json, os, sys
import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ─── orbital parameters ──────────────────────────────────────────────
P_DAYS = 44.41
ECC = 0.007
K1 = 80.355        # km/s
K1_ERR = 0.323

# ─── primary mass (extinction-corrected HR diagram) ──────────────────
M1_BEST = 4.40     # Msun
M1_SIGMA = 1.32    # 30% fractional uncertainty
M1_CAT = 3.16      # catalogue value for comparison

# ─── mass function ───────────────────────────────────────────────────
N_DRAWS = 200_000
NS_MAX = 3.0       # Msun — conservative NS ceiling
BH_THRESHOLD = 5.0 # Msun — upper edge of mass gap


def mass_function(P, K1, e):
    """Spectroscopic mass function f(M) in solar masses."""
    return 1.0385e-7 * (1 - e**2)**1.5 * K1**3 * P


def solve_m2_min(fM, M1):
    """Solve (M2 sin i)^3 / (M1+M2)^2 = f(M) for M2 at i=90°."""
    def eq(m2):
        return m2**3 / (M1 + m2)**2 - fM
    if eq(0.01) > 0:
        return 0.01
    return brentq(eq, 0.01, 500.0)


def main():
    print('=== Bayesian Mass Posterior for HD 338492 ===\n')

    # Mass function
    fM = mass_function(P_DAYS, K1, ECC)
    print(f'  f(M) = {fM:.4f} Msun')
    print(f'  Absolute mass floor (M1→0): {fM:.2f} Msun')

    # M2_min at various M1
    m2_cat = solve_m2_min(fM, M1_CAT)
    m2_best = solve_m2_min(fM, M1_BEST)
    m2_low = solve_m2_min(fM, M1_BEST * 0.7)
    m2_high = solve_m2_min(fM, M1_BEST * 1.3)

    print(f'\n  M2_min (M1={M1_CAT:.2f} cat)  = {m2_cat:.2f} Msun')
    print(f'  M2_min (M1={M1_BEST:.2f} HR)   = {m2_best:.2f} Msun')
    print(f'  M2_min range [{m2_low:.2f}, {m2_high:.2f}] Msun')

    # Monte Carlo inclination debiasing
    print(f'\n  Running MC ({N_DRAWS:,} draws) ...')
    rng = np.random.default_rng(42)

    # Draw M1 from lognormal (positive-definite)
    m1_draws = rng.lognormal(
        mean=np.log(M1_BEST) - 0.5 * (M1_SIGMA / M1_BEST)**2,
        sigma=M1_SIGMA / M1_BEST,
        size=N_DRAWS
    )

    # Draw isotropic inclinations: P(i) ∝ sin(i)
    cos_i = rng.uniform(0, 1, N_DRAWS)
    sin_i = np.sqrt(1 - cos_i**2)

    m2_draws = []
    for m1, si in zip(m1_draws, sin_i):
        # f(M) = (M2 sin i)^3 / (M1+M2)^2
        # Need to solve for M2 given sin_i
        fM_eff = fM / si**3  # effective mass function at this inclination
        try:
            m2 = brentq(lambda m2: m2**3 / (m1 + m2)**2 - fM_eff,
                        0.01, 5000.0)
            m2_draws.append(m2)
        except (ValueError, RuntimeError):
            continue

    m2_draws = np.array(m2_draws)
    n_valid = len(m2_draws)
    print(f'  Valid samples: {n_valid:,}/{N_DRAWS:,}')

    # Statistics
    median = np.median(m2_draws)
    ci68 = np.percentile(m2_draws, [16, 84])
    ci90 = np.percentile(m2_draws, [5, 95])
    p_bh = 100 * np.mean(m2_draws > BH_THRESHOLD)
    p_above_ns = 100 * np.mean(m2_draws > NS_MAX)
    p_mg = 100 * np.mean((m2_draws > NS_MAX) & (m2_draws <= BH_THRESHOLD))
    p_ns = 100 * np.mean(m2_draws <= NS_MAX)

    print(f'\n  Results:')
    print(f'    M2 median     = {median:.2f} Msun')
    print(f'    68% CI        = [{ci68[0]:.2f}, {ci68[1]:.2f}] Msun')
    print(f'    90% CI        = [{ci90[0]:.2f}, {ci90[1]:.2f}] Msun')
    print(f'    P(>NS ceiling) = {p_above_ns:.1f}%')
    print(f'    P(BH > 5)     = {p_bh:.1f}%')
    print(f'    P(mass gap)   = {p_mg:.1f}%')
    print(f'    P(NS)         = {p_ns:.1f}%')

    # Save results
    basedir = os.path.dirname(__file__)
    results = {
        'mass_function_msun': round(fM, 4),
        'M1_best': M1_BEST,
        'M1_sigma': M1_SIGMA,
        'M2_min_catalog_M1': round(m2_cat, 2),
        'M2_min_HR_M1': round(m2_best, 2),
        'M2_min_range': [round(m2_low, 2), round(m2_high, 2)],
        'M2_absolute_floor': round(fM, 2),
        'MC_draws': N_DRAWS,
        'MC_valid': n_valid,
        'M2_median': round(median, 2),
        'M2_68ci': [round(ci68[0], 2), round(ci68[1], 2)],
        'M2_90ci': [round(ci90[0], 2), round(ci90[1], 2)],
        'P_above_NS_percent': round(p_above_ns, 1),
        'P_BH_percent': round(p_bh, 1),
        'P_massgap_percent': round(p_mg, 1),
        'P_NS_percent': round(p_ns, 1),
    }
    outpath = os.path.join(basedir, '..', 'results',
                           'mass_posterior_results.json')
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2)
    print(f'\n  ✓ Saved: {outpath}')

    # Figure
    figpath = os.path.join(basedir, '..', 'paper', 'figures',
                           'fig_mass_posterior.pdf')

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Left: M2_min vs M1
    m1_grid = np.linspace(0.5, 10, 200)
    m2min_grid = [solve_m2_min(fM, m1) for m1 in m1_grid]
    ax1.plot(m1_grid, m2min_grid, 'b-', lw=2)
    ax1.axhline(BH_THRESHOLD, color='red', ls='--', alpha=0.7,
                label=r'$M_2 = 5\,M_\odot$ (BH threshold)')
    ax1.axhline(NS_MAX, color='orange', ls='--', alpha=0.7,
                label=r'$M_2 = 3\,M_\odot$ (NS ceiling)')
    ax1.axhline(fM, color='gray', ls=':', alpha=0.5,
                label=f'$f(M) = {fM:.2f}\\,M_\\odot$ (absolute floor)')
    ax1.axvspan(M1_BEST * 0.7, M1_BEST * 1.3, alpha=0.15, color='blue',
                label=f'$M_1 = {M1_BEST:.1f} \\pm 30\\%$')
    ax1.axvline(M1_BEST, color='blue', ls=':', alpha=0.5)
    ax1.set_xlabel(r'$M_1$ ($M_\odot$)', fontsize=12)
    ax1.set_ylabel(r'$M_{2,\mathrm{min}}$ ($M_\odot$)', fontsize=12)
    ax1.set_title(r'Minimum companion mass vs $M_1$')
    ax1.legend(fontsize=8, loc='upper left')
    ax1.set_xlim(0.5, 10)
    ax1.set_ylim(2, 12)

    # Right: mass posterior histogram
    bins = np.linspace(3, 50, 100)
    ax2.hist(m2_draws[m2_draws < 50], bins=bins, density=True,
             color='steelblue', alpha=0.7, edgecolor='navy', lw=0.3)
    ax2.axvline(BH_THRESHOLD, color='red', ls='--', lw=2,
                label=f'BH threshold (5 $M_\\odot$)')
    ax2.axvline(median, color='black', ls='-', lw=2,
                label=f'Median = {median:.1f} $M_\\odot$')
    ax2.axvspan(ci68[0], ci68[1], alpha=0.15, color='green',
                label=f'68% CI [{ci68[0]:.1f}, {ci68[1]:.1f}]')
    ax2.text(0.97, 0.95, f'$P(\\mathrm{{BH}}) = {p_bh:.1f}\\%$',
             transform=ax2.transAxes, ha='right', va='top', fontsize=14,
             bbox=dict(boxstyle='round', fc='lightyellow'))
    ax2.set_xlabel(r'$M_2$ ($M_\odot$)', fontsize=12)
    ax2.set_ylabel('Probability density', fontsize=12)
    ax2.set_title(f'Inclination-debiased mass posterior ({N_DRAWS//1000}K draws)')
    ax2.legend(fontsize=8)
    ax2.set_xlim(3, 50)

    fig.suptitle('HD 338492 — Mass Constraints', fontsize=13,
                 fontweight='bold')
    fig.tight_layout()
    fig.savefig(figpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  → {figpath}')

    print('\n=== Mass posterior complete ===')


if __name__ == '__main__':
    main()
