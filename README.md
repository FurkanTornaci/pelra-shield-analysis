# Human factors in Indonesian traditional wooden (Pelra) ship accidents: data and code

This repository holds the data and the analysis behind:

> Ardhi, E.; Turan, O.; Kurt, R.; Arif, M.; Tornaci, F. *Human Factors in Traditional Wooden Ship Accidents: An
> Empirical SHIELD Taxonomy Analysis of Indonesian Pelra Vessels.* Applied Sciences (under review).

One script, `analyse.py`, reads the data in `data/` and produces every number, table and data figure in the paper.
Nothing is edited by hand.

## The study in brief

| Question | Data | Method |
|---|---|---|
| RQ1. How does casualty severity in Pelra compare with other fleets? | 215 accident occurrences investigated by Indonesia's National Transportation Safety Committee (NTSC/KNKT), 2007–2022; 208 in the seven commercial fleet categories with at least five occurrences | Kruskal–Wallis, Dunn and Mann–Whitney tests |
| RQ2. Which SHIELD codes and layers recur? | The 23 Pelra cases with a final report: 114 coded instances of 32 SHIELD codes | Counts of cases and instances, weighted by impact rating |
| RQ3. Are SHIELD categories associated with fatalities? | The same 23 cases | Exact Mann–Whitney tests with Benjamini–Hochberg adjustment |
| RQ4. Which interventions should come first? | Ten intervention packages mapped to SHIELD codes, plus support votes from 15 stakeholders | Priority Weight, with six sensitivity analyses |

## Reproduce the results

Requires Python 3.10 or later.

```bash
pip install -r requirements.txt
python analyse.py
```

This takes about a minute and rewrites `outputs/` and `figures/`. The committed copies are the ones in the paper.

## Contents

```
analyse.py         the complete analysis
data/              input data (see data/README.md for every column)
outputs/           tables (CSV) and all statistics (results.json)
figures/           the data figures of the paper (PNG, 300 dpi)
docs/methods.md    the methods in brief, and how the outputs map to the paper
tools/             make_basemap.py, which built the Figure 11 base map from Natural Earth (run once)
```

## Where each result in the paper comes from

| Paper | File |
|---|---|
| Table 1, distribution by layer | `outputs/table1_layers.csv` |
| Table 3, the seven fleet categories, and the Dunn contrasts | `outputs/table3_table10_fleet.csv`, `outputs/dunn_pelra_vs_comparators.csv` |
| Table 4, codes by accident type | `outputs/table4_codes_by_accident_type.csv` |
| Table 5, intervention Priority Weights | `outputs/table5_packages.csv` |
| Table 7, the 215 occurrences | `outputs/S1_accident_occurrences.csv` |
| Table 8, the 30 Pelra cases | `outputs/S2_pelra_30_cases.csv` |
| Table 9, SHIELD codebook | `data/shield_codebook_stroeve2023.csv` |
| Table 10, all 12 fleet categories | `outputs/table3_table10_fleet.csv` |
| Table 11, ranking sensitivity | `outputs/rq4_sensitivity_ranks.csv` |
| Table 12, SHIELD codes: cases, instances and weights | `outputs/table12_codes.csv` |
| Table 13, Priority Weight calculation by code | `outputs/table12_codes.csv` (W_i), `outputs/S4_intervention_audit.csv` |
| Table 14, factor–fatality tests | `outputs/rq3_macro_associations.csv` |
| Supplementary Tables S1–S4 | `outputs/S1_*.csv` to `outputs/S4_*.csv` |
| Figures 1, 2 and 4–18 | `figures/figure<N>_*.png` |
| Every statistic quoted in the text | `outputs/results.json` |

The script does not generate:

- Table 2 (the consultation participants)
- Table 6 (the crosswalk to NTSC recommendations)
- Figure 3 (the SHIELD structure, adapted from Stroeve et al., 2023)

## Data sources

- **Accident occurrences:** compiled from the public investigation reports of the KNKT
  (https://knkt.go.id/investigasi). They form a complete census of the ship accident investigation reports
  published for 2007–2022. Three fatality counts and four location labels were corrected against the final reports;
  the corrections are listed in `data/fatality_corrections.csv` and `data/area_corrections.csv`.
- **SHIELD coding:** coded by the authors from the final reports. Each coded instance carries the supporting
  evidence and a four-point impact rating.
- **Case locations:** the recorded coordinates of the 30 Pelra cases (`data/pelra_case_locations.csv`). Two
  were corrected: the latitude sign for Case 25 and the missing decimal coordinates for Case 27.
- **Base map (Figure 11):** Natural Earth 1:10m countries and lakes (public domain), naturalearthdata.com,
  reduced to the map window by `tools/make_basemap.py`.
- **SHIELD definitions:** quoted from Stroeve et al. (2023), *Safety* 9(1), 14,
  https://doi.org/10.3390/safety9010014 (CC BY 4.0).

## Licence

- **Code** (`analyse.py`): MIT; see `LICENSE`.
- **Data, outputs and figures:** CC BY 4.0; see `LICENSE-DATA.md`.

Please cite the paper when you use them. `CITATION.cff` gives the citation in machine-readable form.

## Contact

Eka Ardhi (corresponding author), Department of Naval Architecture, Ocean and Marine Engineering, University of
Strathclyde: eka.ardhi@strath.ac.uk
