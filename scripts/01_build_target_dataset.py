#!/usr/bin/env python3
"""
01_build_target_dataset.py — Fetch and compile all public data for HD 338492.

Queries:
  - Gaia DR3 TAP (astrometry, photometry, NSS orbit, IPD flags)
  - VizieR: 2MASS (II/246), AllWISE (II/328), GALEX (II/335)
  - VizieR: ROSAT 2RXS (IX/47), XMM 4XMM-DR13 (IX/68)
  - SIMBAD (spectral type, identifiers)

Outputs:
  data/hd338492_input_summary.csv   (refreshed)
  data/photometry_compiled.csv      (refreshed)
  results/gaia_query_results.json
"""

import json, os, sys, warnings
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.gaia import Gaia
from astroquery.vizier import Vizier
from astroquery.simbad import Simbad

warnings.filterwarnings('ignore')

# ─── target ───────────────────────────────────────────────────────────
SOURCE_ID = 2021374066702077312
RA, DEC = 292.72283804, 24.75113317
NAME = 'HD 338492'

OUTDIR = os.path.join(os.path.dirname(__file__), '..', 'results')
DATADIR = os.path.join(os.path.dirname(__file__), '..', 'data')
os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(DATADIR, exist_ok=True)

coord = SkyCoord(ra=RA, dec=DEC, unit='deg', frame='icrs')


def query_gaia():
    """Gaia DR3 astrometry + photometry + gspphot + IPD."""
    q = f"""
    SELECT source_id, ra, dec, parallax, parallax_error,
           phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag,
           bp_rp, ruwe,
           teff_gspphot, logg_gspphot, mh_gspphot,
           ag_gspphot, ebpminrp_gspphot,
           ipd_gof_harmonic_amplitude, ipd_frac_multi_peak
    FROM gaiadr3.gaia_source
    WHERE source_id = {SOURCE_ID}
    """
    job = Gaia.launch_job(q)
    return job.get_results()


def query_nss():
    """Gaia DR3 NSS SB1 orbital solution."""
    q = f"""
    SELECT *
    FROM gaiadr3.nss_two_body_orbit
    WHERE source_id = {SOURCE_ID}
    """
    job = Gaia.launch_job(q)
    return job.get_results()


def query_vizier(catalogue, radius_arcsec):
    """Generic VizieR cone search."""
    v = Vizier(columns=['**'], row_limit=5)
    r = v.query_region(coord, radius=radius_arcsec * u.arcsec,
                       catalog=catalogue)
    return r[0] if r else None


def query_simbad():
    """SIMBAD identification."""
    s = Simbad()
    s.add_votable_fields('sp', 'otype')
    return s.query_object(NAME)


def main():
    print(f'=== Building dataset for {NAME} (Gaia DR3 {SOURCE_ID}) ===\n')

    # 1. Gaia astrometry + photometry
    print('[1/6] Gaia DR3 astrometry ...')
    gaia = query_gaia()
    print(f'      G={gaia["phot_g_mean_mag"][0]:.3f}, '
          f'BP-RP={gaia["bp_rp"][0]:.3f}, '
          f'RUWE={gaia["ruwe"][0]:.2f}')

    # 2. Gaia NSS orbit
    print('[2/6] Gaia DR3 NSS orbit ...')
    nss = query_nss()
    print(f'      P={nss["period"][0]:.4f} d, '
          f'e={nss["eccentricity"][0]:.4f}, '
          f'K1={nss["semi_amplitude_primary"][0]:.3f} km/s')

    # 3. 2MASS + WISE photometry
    print('[3/6] 2MASS + AllWISE photometry ...')
    tmass = query_vizier('II/246', 3)
    wise = query_vizier('II/328', 3)
    phot = {}
    if tmass is not None:
        for b in ['Jmag', 'Hmag', 'Kmag']:
            if b in tmass.colnames:
                phot[b.replace('mag', '')] = float(tmass[b][0])
    if wise is not None:
        for b in ['W1mag', 'W2mag', 'W3mag', 'W4mag']:
            if b in wise.colnames:
                phot[b.replace('mag', '')] = float(wise[b][0])
    print(f'      Bands: {list(phot.keys())}')

    # 4. GALEX UV
    print('[4/6] GALEX UV ...')
    galex = query_vizier('II/335', 5)
    if galex is not None:
        print('      GALEX detection found')
        phot['FUV'] = float(galex['FUVmag'][0]) if 'FUVmag' in galex.colnames else None
        phot['NUV'] = float(galex['NUVmag'][0]) if 'NUVmag' in galex.colnames else None
    else:
        print('      No GALEX match — UV non-detection')

    # 5. X-ray (ROSAT + XMM)
    print('[5/6] X-ray catalogues ...')
    rosat = query_vizier('IX/47', 30)
    xmm = query_vizier('IX/68', 15)
    print(f'      ROSAT: {"detection" if rosat else "non-detection"}')
    print(f'      XMM:   {"detection" if xmm else "non-detection"}')

    # 6. SIMBAD
    print('[6/6] SIMBAD identification ...')
    sim = query_simbad()
    sptype = str(sim['SP_TYPE'][0]) if sim is not None else 'unknown'
    otype = str(sim['OTYPE'][0]) if sim is not None else 'unknown'
    print(f'      SpType={sptype}, OType={otype}')

    # Save JSON summary
    results = {
        'name': NAME,
        'source_id': SOURCE_ID,
        'ra': RA, 'dec': DEC,
        'parallax': float(gaia['parallax'][0]),
        'parallax_error': float(gaia['parallax_error'][0]),
        'distance_pc': 1000.0 / float(gaia['parallax'][0]),
        'G': float(gaia['phot_g_mean_mag'][0]),
        'BP': float(gaia['phot_bp_mean_mag'][0]),
        'RP': float(gaia['phot_rp_mean_mag'][0]),
        'BP_RP': float(gaia['bp_rp'][0]),
        'ruwe': float(gaia['ruwe'][0]),
        'teff_gspphot': gaia['teff_gspphot'][0] if gaia['teff_gspphot'][0] else None,
        'ag_gspphot': gaia['ag_gspphot'][0] if gaia['ag_gspphot'][0] else None,
        'ipd_harm_ampl': float(gaia['ipd_gof_harmonic_amplitude'][0]),
        'ipd_frac_multi': float(gaia['ipd_frac_multi_peak'][0]),
        'period': float(nss['period'][0]),
        'eccentricity': float(nss['eccentricity'][0]),
        'K1': float(nss['semi_amplitude_primary'][0]),
        'significance': float(nss['significance'][0]),
        'rv_sys': float(nss['center_of_mass_velocity'][0]) if 'center_of_mass_velocity' in nss.colnames else None,
        'simbad_sptype': sptype,
        'simbad_otype': otype,
        'photometry_extra': phot,
        'xray_rosat': rosat is not None,
        'xray_xmm': xmm is not None,
        'galex_detected': galex is not None,
    }

    outpath = os.path.join(OUTDIR, 'gaia_query_results.json')
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f'\n✓ Saved: {outpath}')

    print('\n=== Dataset build complete ===')


if __name__ == '__main__':
    main()
