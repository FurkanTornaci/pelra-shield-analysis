# Data

All files are UTF-8. Case numbers (`Case N`, `case_no`, `table8_case`) are the Pelra case numbers used in the paper
(Table 8).

## `ntsc_maritime_accidents_2007_2022.xlsx` and `.csv`

There is one row per accident occurrence investigated by the KNKT (NTSC), 2007–2022 (215 rows). The CSV is an
export of the workbook.

| Column | Meaning |
|---|---|
| Row ID | Row number (1–215, chronological) |
| Year | Year of the accident |
| Case | Pelra case number (`Case 1` … `Case 30`); blank for other fleets |
| Accident Category | Capsizing, Collision, Capsizing after Collision, Fire, Explosion, Grounding or Other |
| Area | Sea, strait, river, lake or port where the accident occurred |
| Class | IMO casualty severity class: VSMC (very serious marine casualty) or SMC (serious marine casualty) |
| Analytical Fleet | Fleet category of the occurrence; in multi-vessel accidents, Pelra if a Pelra vessel was involved |
| Final Fatalities | Persons on board (crew and passengers) who died, plus those missing and presumed dead |
| Fatal Status | Fatal if Final Fatalities > 0, otherwise Non-fatal |
| Consequence | Damage category recorded in the report |
| Factor | Main factor recorded in the report (Human, Technical, Weather; "-" where none is given) |
| Report | KNKT report status. A released report stays "In progress" (open) until the responses of the parties involved have been received or it has been formally approved, then becomes "Final" |

The fleet comparison uses the seven commercial fleet categories with at least five occurrences (208 rows). The five
categories with fewer occurrences (landing craft, barge, fishing vessel, patrol vessel, tug; 7 rows) are kept for
description only.

A vessel is classified as **Pelra** (traditional shipping) when it is built by traditional methods, typically without
formal naval-architecture drawings or classification approval, predominantly of wood, and certified by the local
authority, acting as port authority, to carry general cargo and passengers within its operating port areas. Some
such vessels are built partly of other materials (`pelra_case_attributes.csv`, column `non_wooden_hull`).

## `shield_coding_instances.csv`

There is one row per coded SHIELD instance (114 rows) in the 23 Pelra cases with a final report. A code can occur
more than once in a case.

| Column | Meaning |
|---|---|
| coding_case | Case number in the coding records (1–23, numbered consecutively) |
| table8_case | The same case's Pelra case number in the paper (Table 8) |
| layer, category, code, label | SHIELD layer, category, factor code and factor label |
| impact_rating | Four-point ordinal rating (1 = minor/minimal, 2 = moderate, 3 = major, 4 = critical), assigned by the authors with advice from classification-society experts, the KNKT and the Ministry of Transportation |
| evidence | The coder's summary of the report evidence for the instance |

`coding_case` and `table8_case` differ because the seven excluded cases are interspersed: coding Case 9 is Table 8
Case 10, and so on. Always join on `table8_case`.

## `coded_case_key.csv`

Maps `coding_case` to `table8_case`, with the vessel names that identify each case in the KNKT reports.

## `pelra_case_attributes.csv`

There is one row per Pelra case (30 rows).

| Column | Meaning |
|---|---|
| case_no | Pelra case number |
| vessel_names_cleaned_data | Vessel(s) involved |
| class_column_J | Classification society of each vessel involved (BKI, NK, …, or Non-klas) |
| gt | Gross tonnage of each vessel, where recorded |
| report_status_cleaned_data | Report status in the authors' detailed records |
| classification_group | Non-Classed Only (all vessels non-classed), Mixed Class (a non-classed and a classed vessel) or Classed Only |
| knkt_report | KNKT report number (coded cases) |
| non_wooden_hull | Hull material where the KNKT report describes a non-wooden hull, with the report page |

## `pelra_case_locations.csv`

There is one row per Pelra case (30 rows); Figure 11 plots these.

| Column | Meaning |
|---|---|
| case_no | Pelra case number |
| vessel | Vessel(s) involved |
| location_recorded | Location as recorded in the authors' accident records |
| latitude_recorded, longitude_recorded | Coordinates as recorded (degrees, minutes, seconds; U/T are Indonesian for N/E) |
| latitude, longitude | Decimal degrees used for the map (south and west negative) |
| water_type | Sea, River, Lake or Port |
| note | Where the decimal coordinates differ from the recorded values, the reason |

Case 19 is a port here (Kelapa Lima port on the Maro River, Merauke). The Area column of the accident workbook gives
"River Maro", so Figure 12, which classifies all 215 occurrences by that column, counts it as River.

## `basemap_naturalearth_10m.geojson`

These are the coastlines, country borders and lakes for Figure 11, from Natural Earth 1:10m (public domain). Each
feature lists rings of [longitude, latitude]. `properties.holes` marks rings that are holes, and `properties.kind` is
`country` or `lake`. `tools/make_basemap.py` built this file and can rebuild it (it needs `pyshp`).

## `packages.csv`

The ten intervention packages: code, name, mapped SHIELD codes (separated by `;`) and `n_support`, the number of the
15 consulted stakeholders who identified the package as highly feasible and urgent.

## `shield_codebook_stroeve2023.csv`

Definitions of the 32 SHIELD factors identified in the sample, quoted from Stroeve et al. (2023), *Safety* 9(1), 14,
Appendix A (CC BY 4.0), with the source table and page for each.

## `fatality_corrections.csv`, `area_corrections.csv`

These are corrections made against the KNKT final reports and already applied in the accident workbook. Each row
gives the original value, the corrected value and the report and page used. `analyse.py` checks that they are in
place.
