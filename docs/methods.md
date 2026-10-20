# Methods in brief

The paper gives the full description; this note explains what `analyse.py` computes.

## Fleet comparison (RQ1)

- **Groups:** the seven commercial fleet categories with at least five occurrences (n = 208). Pelra (n = 30) is compared with the six other categories (n = 178).
- **Tests:**
  - a Kruskal–Wallis test across the seven categories (effect size η²H = (H − k + 1)/(n − k))
  - Dunn's z-tests for the six Pelra-versus-comparator contrasts (tie-corrected, Bonferroni ×6)
  - a Mann–Whitney U test of Pelra against the six comparators (normal approximation with tie and continuity correction; r = |z|/√N; Hodges–Lehmann shift with a distribution-free 95% CI)
  - Fisher's exact test for the proportion of fatal occurrences
- **Summaries:** proportions carry Wilson 95% CIs, and quartiles use linear interpolation.

## SHIELD frequencies (RQ2)

- **Cases:** the number of the 23 cases in which a code is present.
- **Instances:** the number of times the code was recorded.
- **Weight W_i:** the sum of the 1–4 impact ratings of the code's instances.

Layer and category totals are built from the same 23 × 32 matrices (`outputs/S3_*`).

## Factor–fatality associations (RQ3)

- **Categories tested:**
  - Operations Planning (LO1, LO2, LO6)
  - Culture (OC1)
  - Resources (OR3, OR5)
  - Safety Management (OS2)
- **Test:** a two-sided Mann–Whitney U test with exact permutation p-values, which remain valid with ties and small groups. Benjamini–Hochberg adjustment is applied over the four tests (Q = 0.10).
- **Sensitivity check:** 95th-percentile winsorisation.

## Intervention Priority Weights (RQ4)

PW_j = S_j × Σ_{i ∈ C_j} F_i · I_i, where:

- **F_i:** the instances of code i.
- **I_i:** the mean impact rating of code i, so that F_i · I_i = W_i.
- **S_j:** 1 + n_support/15.

Six sensitivity analyses re-rank the packages:

1. case counts with a median tier per code
2. frequency only
3. geometric weights (1, 2, 4, 8)
4. the mean instead of the sum over mapped codes
5. case coverage
6. no S_j

## Sensitivity to vessel identification

The main analyses are repeated without the three coded vessels that the KNKT reports describe as non-wooden (Cases 7,
26 and 28). The results are in `results.json` under `sensitivity_non_wooden_excluded`.

## Map (Figure 11)

Each Pelra case is plotted at its recorded coordinates, coloured by water type. Markers that would overlap a
neighbour are offset and joined to the true location by a line. The base map uses an equirectangular (plate carrée)
projection, which near the equator distorts little.
