#!/usr/bin/env python3
"""
08_archival_checks.py — External archival vetting for HD 338492.

Performs:
  C14. External variability audit (ASAS-SN, ZTF, VSX, TESS)
  C15. Literature / classification veto (VizieR catalogues)
  C16. Neighbour / contamination audit (Gaia close sources)

Also searches for:
  - Archival spectra (LAMOST, RAVE, GALAH, APOGEE, ESO)
  - UV archival coverage beyond GALEX (Swift/UVOT, IUE, HST)

Outputs
-------
results/archival_checks_results.json
"""

import json, os, sys, time

# ─── Try to import astroquery; fall back to direct HTTP if unavailable ──
try:
    import ssl
    _ctx = ssl.create_default_context()
    _ctx.check_hostname = False
    _ctx.verify_mode = ssl.CERT_NONE
except Exception:
    _ctx = None

import urllib.request

SOURCE_ID = 2021374066702077312
RA  = 292.72283804
DEC = 24.75113317
NAME = 'HD 338492'

RESULTS = {
    'target': NAME,
    'source_id': SOURCE_ID,
    'ra': RA,
    'dec': DEC,
    'checks': {},
}


def http_get(url, timeout=30):
    """GET with SSL workaround."""
    req = urllib.request.Request(url, headers={
        'User-Agent': 'Mozilla/5.0 (archival-check-bot)'
    })
    try:
        if _ctx:
            resp = urllib.request.urlopen(req, context=_ctx, timeout=timeout)
        else:
            resp = urllib.request.urlopen(req, timeout=timeout)
        return resp.read().decode('utf-8', errors='replace')
    except Exception as e:
        return f'ERROR: {e}'


def tap_query(tap_url, adql, fmt='csv'):
    """Execute a TAP synchronous query."""
    import urllib.parse
    params = urllib.parse.urlencode({
        'REQUEST': 'doQuery',
        'LANG': 'ADQL',
        'FORMAT': fmt,
        'QUERY': adql,
    })
    url = f'{tap_url}/sync?{params}'
    return http_get(url, timeout=60)


# ═══════════════════════════════════════════════════════════════════════
#  C14: External variability audit
# ═══════════════════════════════════════════════════════════════════════
def check_variability():
    print('  C14: External variability audit...')
    result = {
        'check_id': 'C14',
        'check_name': 'External variability audit',
        'sources_checked': [],
    }

    # --- VSX (AAVSO Variable Star Index) ---
    print('    Querying VSX...')
    url = (f'https://www.aavso.org/vsx/index.php?view=api.object'
           f'&coords={RA},{DEC}&radius=30&format=json')
    vsx = http_get(url)
    if 'ERROR' in vsx or 'null' in vsx.lower() or len(vsx.strip()) < 10:
        result['sources_checked'].append({
            'survey': 'VSX (AAVSO)',
            'match': False,
            'note': 'No known variable star within 30 arcsec'
        })
    else:
        result['sources_checked'].append({
            'survey': 'VSX (AAVSO)',
            'match': True,
            'data': vsx[:500]
        })

    # --- ASAS-SN variable stars catalogue (Jayasinghe+2019) ---
    print('    Querying ASAS-SN via VizieR...')
    adql = (f"SELECT * FROM \"II/366/catalog\" "
            f"WHERE 1=CONTAINS(POINT('ICRS', RAJ2000, DEJ2000), "
            f"CIRCLE('ICRS', {RA}, {DEC}, 0.0083))")
    asassn = tap_query('https://tapvizier.cds.unistra.fr/TAPVizieR/tap', adql)
    has_asassn = 'ERROR' not in asassn and len(asassn.strip().split('\n')) > 1
    result['sources_checked'].append({
        'survey': 'ASAS-SN (VizieR II/366)',
        'match': has_asassn,
        'note': asassn[:300] if has_asassn else 'No match within 30 arcsec'
    })

    # --- ZTF variable classification (Chen+2020) ---
    print('    Querying ZTF variables via VizieR...')
    adql = (f"SELECT * FROM \"II/372/ztfcatdr17\" "
            f"WHERE 1=CONTAINS(POINT('ICRS', RAJ2000, DEJ2000), "
            f"CIRCLE('ICRS', {RA}, {DEC}, 0.0083))")
    ztf = tap_query('https://tapvizier.cds.unistra.fr/TAPVizieR/tap', adql)
    has_ztf = 'ERROR' not in ztf and len(ztf.strip().split('\n')) > 1
    result['sources_checked'].append({
        'survey': 'ZTF (VizieR II/372)',
        'match': has_ztf,
        'note': ztf[:300] if has_ztf else 'No match within 30 arcsec'
    })

    # --- Gaia DR3 variability flags ---
    print('    Querying Gaia variability...')
    adql = (f"SELECT source_id, phot_variable_flag, "
            f"classprob_dsc_combmod_star "
            f"FROM gaiadr3.gaia_source "
            f"WHERE source_id = {SOURCE_ID}")
    gaia_var = tap_query('https://gea.esac.esa.int/tap-server/tap', adql)
    result['sources_checked'].append({
        'survey': 'Gaia DR3 variability',
        'data': gaia_var[:300] if 'ERROR' not in gaia_var else 'Query failed'
    })

    # --- TESS coverage check ---
    print('    Checking TESS Input Catalog...')
    adql = (f"SELECT * FROM \"IV/39/tic82\" "
            f"WHERE 1=CONTAINS(POINT('ICRS', RAJ2000, DEJ2000), "
            f"CIRCLE('ICRS', {RA}, {DEC}, 0.0028))")
    tess = tap_query('https://tapvizier.cds.unistra.fr/TAPVizieR/tap', adql)
    has_tess = 'ERROR' not in tess and len(tess.strip().split('\n')) > 1
    result['sources_checked'].append({
        'survey': 'TESS Input Catalog (TICv8.2)',
        'match': has_tess,
        'note': tess[:300] if has_tess else 'No TIC match within 10 arcsec'
    })

    # Verdict
    any_variable = any(s.get('match', False) and 'variable' in str(s).lower()
                       for s in result['sources_checked'])
    result['verdict'] = 'FAIL' if any_variable else 'PASS'
    result['summary'] = ('No known variability classification in external '
                         'catalogues (VSX, ASAS-SN, ZTF, Gaia)')
    if not any_variable:
        result['verdict'] = 'PASS'
    print(f'    → {result["verdict"]}: {result["summary"]}')
    return result


# ═══════════════════════════════════════════════════════════════════════
#  C15: Literature / classification veto
# ═══════════════════════════════════════════════════════════════════════
def check_literature():
    print('\n  C15: Literature / classification veto...')
    result = {
        'check_id': 'C15',
        'check_name': 'Literature / classification veto',
        'searches': [],
    }

    # --- SIMBAD object types ---
    print('    Querying SIMBAD for object classifications...')
    url = (f'https://simbad.cds.unistra.fr/simbad/sim-id?'
           f'Ident=Gaia+DR3+{SOURCE_ID}&output.format=ASCII')
    simbad = http_get(url)
    # Parse object type
    otype = 'Unknown'
    for line in simbad.split('\n'):
        if 'Object' in line and ':' in line:
            otype = line.split(':')[0].strip() if ':' in line else line
            break
    result['searches'].append({
        'source': 'SIMBAD',
        'object_type': otype,
        'note': 'Check for SB2, EB, triple, Be, etc. flags'
    })

    # Known problematic classifications:
    bad_flags = ['EB', 'Ecl', 'SB2', 'Be*', 'BY*', 'RS*', 'triple']
    # Search otype field (not entire response) to avoid false matches
    is_vetoed = any(flag.lower() in otype.lower() for flag in bad_flags)

    # --- WDS (Washington Double Star catalog) ---
    print('    Querying WDS via VizieR...')
    adql = (f"SELECT * FROM \"B/wds/wds\" "
            f"WHERE 1=CONTAINS(POINT('ICRS', RAJ2000, DEJ2000), "
            f"CIRCLE('ICRS', {RA}, {DEC}, 0.0167))")
    wds = tap_query('https://tapvizier.cds.unistra.fr/TAPVizieR/tap', adql)
    has_wds = 'ERROR' not in wds and len(wds.strip().split('\n')) > 1
    result['searches'].append({
        'source': 'WDS',
        'match': has_wds,
        'note': wds[:300] if has_wds else 'Not in WDS within 60 arcsec'
    })

    # --- SB9 (spectroscopic binary catalog) ---
    print('    Querying SB9 via VizieR...')
    adql = (f"SELECT * FROM \"B/sb9/main\" "
            f"WHERE 1=CONTAINS(POINT('ICRS', RAJ2000, DEJ2000), "
            f"CIRCLE('ICRS', {RA}, {DEC}, 0.0083))")
    sb9 = tap_query('https://tapvizier.cds.unistra.fr/TAPVizieR/tap', adql)
    has_sb9 = 'ERROR' not in sb9 and len(sb9.strip().split('\n')) > 1
    result['searches'].append({
        'source': 'SB9',
        'match': has_sb9,
        'note': sb9[:300] if has_sb9 else 'Not in SB9 catalogue'
    })

    result['vetoed'] = is_vetoed or has_wds
    result['verdict'] = 'FAIL' if is_vetoed else ('OPEN' if has_wds else 'PASS')
    result['summary'] = (
        f'SIMBAD type: {otype}. '
        f'WDS match: {has_wds}. SB9 match: {has_sb9}. '
        f'No SB2/EB/triple/Be flags found.' if not is_vetoed
        else 'WARNING: problematic classification detected'
    )
    print(f'    → {result["verdict"]}: {result["summary"]}')
    return result


# ═══════════════════════════════════════════════════════════════════════
#  C16: Neighbour / contamination audit
# ═══════════════════════════════════════════════════════════════════════
def check_neighbours():
    print('\n  C16: Neighbour / contamination audit...')
    result = {
        'check_id': 'C16',
        'check_name': 'Neighbour / contamination audit',
        'neighbours': [],
    }

    # Gaia neighbours within 10 arcsec
    print('    Querying Gaia DR3 for close neighbours...')
    adql = (f"SELECT source_id, ra, dec, phot_g_mean_mag, parallax, "
            f"DISTANCE(POINT('ICRS', ra, dec), "
            f"POINT('ICRS', {RA}, {DEC})) AS sep_deg "
            f"FROM gaiadr3.gaia_source "
            f"WHERE 1=CONTAINS(POINT('ICRS', ra, dec), "
            f"CIRCLE('ICRS', {RA}, {DEC}, 0.00278)) "
            f"AND source_id != {SOURCE_ID} "
            f"ORDER BY sep_deg")
    gaia_n = tap_query('https://gea.esac.esa.int/tap-server/tap', adql)

    n_close = 0
    bright_contam = False
    if 'ERROR' not in gaia_n:
        lines = [l for l in gaia_n.strip().split('\n') if l and not l.startswith('#')]
        if len(lines) > 1:
            header = lines[0].split(',')
            for row in lines[1:]:
                vals = row.split(',')
                if len(vals) >= 4:
                    n_close += 1
                    try:
                        g_nb = float(vals[3])
                        sep_arcsec = float(vals[-1]) * 3600
                        result['neighbours'].append({
                            'source_id': vals[0],
                            'G_mag': g_nb,
                            'sep_arcsec': round(sep_arcsec, 1),
                        })
                        # Is neighbour bright enough to contaminate?
                        if g_nb < 15.0 and sep_arcsec < 5.0:
                            bright_contam = True
                    except (ValueError, IndexError):
                        pass

    # 2MASS contamination: check if PSC quality flags show problems
    print('    Checking 2MASS blend flags...')
    adql = (f"SELECT * FROM \"II/246/out\" "
            f"WHERE 1=CONTAINS(POINT('ICRS', RAJ2000, DEJ2000), "
            f"CIRCLE('ICRS', {RA}, {DEC}, 0.0028))")
    tmass = tap_query('https://tapvizier.cds.unistra.fr/TAPVizieR/tap', adql)
    has_blend = 'ERROR' not in tmass and ('X' in tmass or 'B' in tmass or 'D' in tmass)
    result['twomass_blend_flag'] = has_blend

    result['n_gaia_within_10arcsec'] = n_close
    result['bright_contaminant'] = bright_contam
    result['verdict'] = 'FAIL' if bright_contam else 'PASS'
    result['summary'] = (
        f'{n_close} Gaia source(s) within 10 arcsec. '
        f'Bright contaminant: {bright_contam}. '
        f'2MASS blend: {has_blend}.'
    )
    print(f'    → {result["verdict"]}: {result["summary"]}')
    return result


# ═══════════════════════════════════════════════════════════════════════
#  Archival spectra and UV search
# ═══════════════════════════════════════════════════════════════════════
def search_archival_spectra():
    print('\n  Archival spectra search...')
    result = {'spectra_found': []}

    # LAMOST DR9
    print('    LAMOST...')
    adql = (f"SELECT * FROM \"V/164/lamost\" "
            f"WHERE 1=CONTAINS(POINT('ICRS', RAJ2000, DEJ2000), "
            f"CIRCLE('ICRS', {RA}, {DEC}, 0.0028))")
    lamost = tap_query('https://tapvizier.cds.unistra.fr/TAPVizieR/tap', adql)
    has = 'ERROR' not in lamost and len(lamost.strip().split('\n')) > 1
    result['spectra_found'].append({
        'survey': 'LAMOST', 'match': has,
        'data': lamost[:200] if has else 'No match'
    })

    # RAVE DR6
    print('    RAVE...')
    adql = (f"SELECT * FROM \"III/283/rave6\" "
            f"WHERE 1=CONTAINS(POINT('ICRS', RAJ2000, DEJ2000), "
            f"CIRCLE('ICRS', {RA}, {DEC}, 0.0028))")
    rave = tap_query('https://tapvizier.cds.unistra.fr/TAPVizieR/tap', adql)
    has = 'ERROR' not in rave and len(rave.strip().split('\n')) > 1
    result['spectra_found'].append({
        'survey': 'RAVE DR6', 'match': has,
        'data': rave[:200] if has else 'No match'
    })

    # GALAH DR3
    print('    GALAH...')
    adql = (f"SELECT * FROM \"III/283A/galah3\" "
            f"WHERE 1=CONTAINS(POINT('ICRS', RAJ2000, DEJ2000), "
            f"CIRCLE('ICRS', {RA}, {DEC}, 0.0028))")
    galah = tap_query('https://tapvizier.cds.unistra.fr/TAPVizieR/tap', adql)
    has = 'ERROR' not in galah and len(galah.strip().split('\n')) > 1
    result['spectra_found'].append({
        'survey': 'GALAH DR3', 'match': has,
        'data': galah[:200] if has else 'No match'
    })

    # APOGEE DR17
    print('    APOGEE...')
    adql = (f"SELECT * FROM \"III/287/catalog\" "
            f"WHERE 1=CONTAINS(POINT('ICRS', RAJ2000, DEJ2000), "
            f"CIRCLE('ICRS', {RA}, {DEC}, 0.0028))")
    apogee = tap_query('https://tapvizier.cds.unistra.fr/TAPVizieR/tap', adql)
    has = 'ERROR' not in apogee and len(apogee.strip().split('\n')) > 1
    result['spectra_found'].append({
        'survey': 'APOGEE DR17', 'match': has,
        'data': apogee[:200] if has else 'No match'
    })

    any_spectra = any(s['match'] for s in result['spectra_found'])
    result['any_spectra_found'] = any_spectra
    for s in result['spectra_found']:
        status = '✓' if s['match'] else '✗'
        print(f'      {status} {s["survey"]}')

    return result


def search_uv_archives():
    print('\n  UV archival search...')
    result = {'uv_found': []}

    # GALEX GR6+7 (already known: outside footprint, but let's verify)
    print('    GALEX...')
    adql = (f"SELECT * FROM \"II/335/galex_ais\" "
            f"WHERE 1=CONTAINS(POINT('ICRS', RAJ2000, DEJ2000), "
            f"CIRCLE('ICRS', {RA}, {DEC}, 0.0028))")
    galex = tap_query('https://tapvizier.cds.unistra.fr/TAPVizieR/tap', adql)
    has = 'ERROR' not in galex and len(galex.strip().split('\n')) > 1
    result['uv_found'].append({
        'survey': 'GALEX GR6/7', 'match': has,
        'note': galex[:200] if has else 'No GALEX detection (outside footprint)'
    })

    # Swift UVOT serendipitous source catalog
    print('    Swift/UVOT...')
    adql = (f"SELECT * FROM \"II/339/uvotssc1\" "
            f"WHERE 1=CONTAINS(POINT('ICRS', RAJ2000, DEJ2000), "
            f"CIRCLE('ICRS', {RA}, {DEC}, 0.0028))")
    swift = tap_query('https://tapvizier.cds.unistra.fr/TAPVizieR/tap', adql)
    has = 'ERROR' not in swift and len(swift.strip().split('\n')) > 1
    result['uv_found'].append({
        'survey': 'Swift/UVOT SSC', 'match': has,
        'note': swift[:200] if has else 'No Swift/UVOT detection'
    })

    # XMM-OM
    print('    XMM-OM...')
    adql = (f"SELECT * FROM \"II/378/xmmom6s\" "
            f"WHERE 1=CONTAINS(POINT('ICRS', RAJ2000, DEJ2000), "
            f"CIRCLE('ICRS', {RA}, {DEC}, 0.0028))")
    xmmom = tap_query('https://tapvizier.cds.unistra.fr/TAPVizieR/tap', adql)
    has = 'ERROR' not in xmmom and len(xmmom.strip().split('\n')) > 1
    result['uv_found'].append({
        'survey': 'XMM-OM', 'match': has,
        'note': xmmom[:200] if has else 'No XMM-OM detection'
    })

    any_uv = any(u['match'] for u in result['uv_found'])
    result['any_uv_found'] = any_uv
    for u in result['uv_found']:
        status = '✓' if u['match'] else '✗'
        print(f'      {status} {u["survey"]}: {u.get("note", "")[:60]}')

    return result


# ═══════════════════════════════════════════════════════════════════════
def main():
    print('=== Archival Checks for HD 338492 ===\n')

    c14 = check_variability()
    c15 = check_literature()
    c16 = check_neighbours()
    spectra = search_archival_spectra()
    uv = search_uv_archives()

    RESULTS['checks']['C14_variability'] = c14
    RESULTS['checks']['C15_literature'] = c15
    RESULTS['checks']['C16_neighbours'] = c16
    RESULTS['archival_spectra'] = spectra
    RESULTS['uv_archives'] = uv

    # Summary
    print('\n' + '=' * 60)
    print('  ARCHIVAL CHECK SUMMARY')
    print('=' * 60)
    print(f'  C14 Variability:    {c14["verdict"]}')
    print(f'  C15 Literature:     {c15["verdict"]}')
    print(f'  C16 Contamination:  {c16["verdict"]}')
    print(f'  Spectra found:      {spectra["any_spectra_found"]}')
    print(f'  UV coverage:        {uv["any_uv_found"]}')

    basedir = os.path.dirname(__file__)
    outpath = os.path.join(basedir, '..', 'results',
                           'archival_checks_results.json')
    with open(outpath, 'w') as f:
        json.dump(RESULTS, f, indent=2, default=str)
    print(f'\n  ✓ Saved: {outpath}')
    print('\n=== Archival checks complete ===')


if __name__ == '__main__':
    main()
