# HD 338492 — Data Notes

## Source Identification
- **HD 338492** = Gaia DR3 2021374066702077312
- SIMBAD object type: SB* (spectroscopic binary)
- SIMBAD spectral type: B9
- Galactic coordinates: (l, b) = (59.30°, 3.01°) — deep in the Galactic plane

## Gaia DR3 NSS SB1 Solution
The orbital solution comes from the Gaia DR3 Non-Single Star (NSS)
catalogue, table `gaiadr3.nss_two_body_orbit`. Solution type: SB1
(single-lined spectroscopic binary).

- The `gspphot` atmospheric parameters (Teff, logg, [M/H], A_G,
  E(BP−RP)) are **not available** for this source — all return NULL.
  This is common for early-type or heavily reddened stars where the
  GSP-Phot pipeline struggles.

## Photometry Sources
- **Gaia DR3**: G, BP, RP from `gaiadr3.gaia_source`
- **2MASS** (catalogue II/246): J, H, Ks from VizieR
- **AllWISE** (catalogue II/328): W1, W2, W3, W4 from VizieR
- **GALEX**: No match within 5 arcsec — no UV detection

## Critical Finding: Reddening
The SIMBAD spectral type B9 implies Teff ≈ 10,500 K, but the
observed BP−RP = 1.039 corresponds to a photometric Teff of only
~5,400 K. This 5,100 K discrepancy is fully explained by interstellar
reddening along the Galactic plane sightline:
- E(BP−RP) ≈ 1.189 mag
- A_G ≈ 2.247 mag
- A_V ≈ 2.85 mag

## No Existing Spectra
This source has no entries in GALAH DR3, RAVE DR6, LAMOST, or
APOGEE. The only radial velocity measurement is from Gaia DR3 NSS.
