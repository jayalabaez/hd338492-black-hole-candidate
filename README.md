# HD 338492 — Compelling Black Hole Candidate from Gaia DR3

> **Joel Ayala-Baez** | Independent Researcher  
> Paper status: Pre-submission draft

## Summary

HD 338492 (Gaia DR3 2021374066702077312) is a single-lined spectroscopic binary (SB1) harbouring a **compelling dormant stellar-mass black hole candidate**. This repository contains all data, analysis scripts, and the manuscript for the associated publication.

### Key Results

| Property | Value |
|---|---|
| Spectral type | B9 (SIMBAD) |
| Period | 44.41 d |
| Eccentricity | 0.007 |
| K₁ | 80.4 km/s |
| f(M) | 2.39 M☉ (exceeds Chandrasekhar limit) |
| M₁ (primary) | 4.40 ± 1.32 M☉ |
| M₂,min (companion) | 6.62 M☉ |
| P(BH) | 98.5% |
| Companion light | Excluded (419× threshold) |
| X-ray | Non-detection (dormant) |
| Alternative scenarios | 4/6 excluded, 2 constrained |

## Repository Structure

```
hd338492-black-hole-candidate/
├── data/                          # Input data and metadata
│   ├── hd338492_input_summary.csv
│   ├── photometry_compiled.csv
│   └── metadata_notes.md
├── scripts/                       # Reproducibility scripts
│   ├── 01_build_target_dataset.py   # Query Gaia, VizieR, SIMBAD
│   ├── 02_fit_sed_extinction.py     # SED analysis + extinction
│   ├── 03_compute_mass_posterior.py  # Mass function + MC posterior
│   ├── 04_companion_exclusion.py    # Luminous companion test
│   ├── 05_alternative_scenarios.py  # Non-BH scenario assessment
│   └── 06_make_figures.py           # Publication figures
├── paper/                         # MNRAS manuscript
│   ├── manuscript.tex
│   ├── references.bib
│   ├── figures/                   # Generated figures (PDF)
│   └── tables/
├── results/                       # Script outputs (JSON)
├── docs/                          # Additional documentation
├── requirements.txt
├── CITATION.cff
└── LICENSE
```

## Quick Start

```bash
# Clone
git clone https://github.com/jayalabaez/hd338492-black-hole-candidate.git
cd hd338492-black-hole-candidate

# Create environment
python -m venv .venv
.venv\Scripts\activate        # Windows
# source .venv/bin/activate   # Linux/Mac

# Install dependencies
pip install -r requirements.txt

# Run analysis (scripts 02–06 are offline; script 01 queries online archives)
python scripts/01_build_target_dataset.py
python scripts/02_fit_sed_extinction.py
python scripts/03_compute_mass_posterior.py
python scripts/04_companion_exclusion.py
python scripts/05_alternative_scenarios.py
python scripts/06_make_figures.py
```

## Data Sources

- **Gaia DR3**: Astrometry, photometry, NSS orbital solution ([Gaia Archive](https://gea.esac.esa.int/archive/))
- **2MASS**: J, H, Ks photometry (via VizieR)
- **AllWISE**: W1–W4 photometry (via VizieR)
- **ROSAT / XMM-Newton**: X-ray upper limits (via VizieR)
- **SIMBAD**: Spectral classification (B9)

## Citation

If you use this work, please cite:

```
Ayala-Baez, Joel (2025). HD 338492 as a Compelling Black Hole Candidate
from Gaia DR3 Astrometric–Spectroscopic Binary Analysis.
```

See [CITATION.cff](CITATION.cff) for machine-readable citation metadata.

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE).
