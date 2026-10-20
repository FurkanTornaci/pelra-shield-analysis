"""Reproduce every quantitative result, table and data figure of the Pelra/SHIELD paper.

Inputs (data/):
  ntsc_maritime_accidents_2007_2022.xlsx
                               215 NTSC-investigated maritime accident occurrences, 2007-2022
  fatality_corrections.csv     three counts corrected against the KNKT reports (already applied in the workbook)
  area_corrections.csv         four area labels corrected against the KNKT reports (already applied)
  shield_coding_instances.csv  114 SHIELD-coded instances in the 23 analysed Pelra cases
                               (with coded report evidence and 1-4 impact ratings)
  coded_case_key.csv           coding-sheet case number (1-23) -> Table 8 case number, vessel names
  packages.csv                 the ten intervention packages, mapped SHIELD codes and stakeholder votes
  pelra_case_attributes.csv    per Pelra case: recorded vessel classes and the derived
                               classification group, KNKT report number, non-wooden hull flag
Outputs: outputs/*.csv, outputs/results.json, figures/*.png
Run:  python analyse.py      (Python 3.10; numpy, pandas, scipy, statsmodels, openpyxl, matplotlib)
"""
import json
import math
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.proportion import proportion_confint

HERE = Path(__file__).resolve().parent
DATA, OUT, FIG = HERE / "data", HERE / "outputs", HERE / "figures"
OUT.mkdir(exist_ok=True)
FIG.mkdir(exist_ok=True)

# Fleet categories with fewer occurrences (barge, tug, patrol vessel, fishing vessel, landing craft) are kept in
# the descriptive inventory only: the comparison focuses on commercial cargo and passenger fleets, and estimates
# from 1-3 occurrences are not stable enough for between-group tests.
MIN_FLEET_N = 5
FLEET_ORDER = ["Pelra", "Ro-Ro Ferry", "General Cargo", "Tanker", "Passenger", "Container", "Bulk Carrier"]
LAYERS = ["Acts", "Preconditions", "Operational Leadership", "Organisation"]
# Category abbreviations as in Figure 3 (Stroeve et al., 2023)
CATEGORY_ABBREV = {
    "Planning And Decision Making": "PD", "Perception": "Pn", "Communicating": "Cn", "Intentional Deviation": "ID",
    "Equipment And Workplace": "EW", "Competence, Skills And Capability": "CS", "Awareness": "As",
    "Physical Environment": "PE", "Physiological Condition": "PC", "Task Leadership": "TL",
    "Operations Planning": "OP", "Personnel Leadership": "PL", "Resources": "Rs", "Culture": "Ce",
    "Economy And Business": "EB", "Safety Management": "SM",
}
CATEGORY_NAME = {
    "Planning And Decision Making": "Planning and Decision Making", "Communicating": "Communication",
    "Equipment And Workplace": "Equipment and Workplace", "Competence, Skills And Capability":
    "Competence, Skills and Capability", "Economy And Business": "Economy and Business",
}
# Macro-categories tested against fatalities (RQ3), defined by SHIELD code
MACRO = {
    "Operations Planning": ["LO1", "LO2", "LO6"],
    "Culture": ["OC1"],
    "Resources": ["OR3", "OR5"],
    "Safety Management": ["OS2"],
}
SEED = 20261020
results = {}


# ---------------------------------------------------------------- statistics helpers
def quartiles(x):
    """Q1 and Q3 by linear interpolation (Hyndman-Fan type 7, the numpy/R/Excel default)."""
    return float(np.percentile(x, 25)), float(np.percentile(x, 75))


def describe(x):
    x = np.asarray(x, float)
    q1, q3 = quartiles(x)
    return dict(n=len(x), mean=float(x.mean()), median=float(np.median(x)), q1=q1, q3=q3, max=float(x.max()),
                total=float(x.sum()), fatal=int((x > 0).sum()), pct_fatal=float(100 * (x > 0).mean()))


def tie_term(values):
    _, t = np.unique(values, return_counts=True)
    return float((t ** 3 - t).sum())


def exact_mwu_p(x, y):
    """Exact two-sided permutation p for the Mann-Whitney U (all C(N, n1) relabellings; handles ties)."""
    allv = np.concatenate([x, y])
    ranks2 = np.round(2 * stats.rankdata(allv)).astype(int)  # doubled mid-ranks are integers
    n1, N = len(x), len(allv)
    smax = int(ranks2.sum())
    dp = np.zeros((n1 + 1, smax + 1))  # dp[k, s]: number of k-subsets with doubled rank sum s
    dp[0, 0] = 1
    for r in ranks2:
        dp[1:, r:] += dp[:-1, :-r].copy()
    sums = np.arange(smax + 1)
    u = sums / 2 - n1 * (n1 + 1) / 2
    mu = n1 * (N - n1) / 2
    u_obs = ranks2[:n1].sum() / 2 - n1 * (n1 + 1) / 2
    extreme = np.abs(u - mu) >= abs(u_obs - mu) - 1e-9
    return float(dp[n1, extreme].sum() / math.comb(N, n1))


def hodges_lehmann(x, y, alpha=0.05):
    """Hodges-Lehmann shift (x - y) with the distribution-free CI based on the normal approximation to U."""
    d = np.sort(np.subtract.outer(np.asarray(x, float), np.asarray(y, float)).ravel())
    n1, n2 = len(x), len(y)
    z = stats.norm.ppf(1 - alpha / 2)
    c = int(math.floor(n1 * n2 / 2 - z * math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)))
    c = max(c, 1)
    return float(np.median(d)), float(d[c - 1]), float(d[n1 * n2 - c])


def mwu(x, y, exact=False):
    """Two-sided Mann-Whitney U, group x vs y.

    Asymptotic p: normal approximation with tie correction and continuity correction (scipy).
    z and r = |z|/sqrt(N) use the tie-corrected variance without continuity correction.
    """
    x, y = np.asarray(x, float), np.asarray(y, float)
    n1, n2 = len(x), len(y)
    N = n1 + n2
    res = stats.mannwhitneyu(x, y, alternative="two-sided", method="asymptotic", use_continuity=True)
    u_x = float(res.statistic)
    sd = math.sqrt(n1 * n2 / 12 * ((N + 1) - tie_term(np.concatenate([x, y])) / (N * (N - 1))))
    z = (u_x - n1 * n2 / 2) / sd
    hl, lo, hi = hodges_lehmann(x, y)
    out = dict(n1=n1, n2=n2, U_x=u_x, U_y=n1 * n2 - u_x, U=min(u_x, n1 * n2 - u_x), z=float(z),
               p_asymptotic=float(res.pvalue), r=abs(z) / math.sqrt(N),
               x=describe(x), y=describe(y), hl_shift=hl, hl_ci95=[lo, hi])
    if exact:
        out["p_exact"] = exact_mwu_p(x, y)
    return out


def kruskal(groups, names):
    groups = [np.asarray(g, float) for g in groups]
    H, p = stats.kruskal(*groups)
    n, k = sum(map(len, groups)), len(groups)
    return dict(H=float(H), df=k - 1, p=float(p), n=n, k=k, eta2_H=(H - k + 1) / (n - k), epsilon2=H / (n - 1),
                groups={nm: describe(g) for nm, g in zip(names, groups)})


def dunn(groups, names, ref, n_contrasts):
    """Dunn's z-tests of `ref` vs every other group (tie-corrected), Bonferroni-adjusted for n_contrasts."""
    allv = np.concatenate(groups)
    ranks = stats.rankdata(allv)
    N = len(allv)
    sizes = [len(g) for g in groups]
    idx = np.cumsum([0] + sizes)
    mean_rank = {nm: float(ranks[idx[i]:idx[i + 1]].mean()) for i, nm in enumerate(names)}
    var0 = N * (N + 1) / 12 - tie_term(allv) / (12 * (N - 1))
    n_of = dict(zip(names, sizes))
    rows = []
    for nm in names:
        if nm == ref:
            continue
        se = math.sqrt(var0 * (1 / n_of[ref] + 1 / n_of[nm]))
        z = (mean_rank[ref] - mean_rank[nm]) / se
        p = 2 * stats.norm.sf(abs(z))
        rows.append(dict(contrast=f"{ref} vs {nm}", mean_rank_ref=mean_rank[ref], mean_rank_other=mean_rank[nm],
                         z=z, p_raw=p, p_bonferroni=min(1.0, p * n_contrasts),
                         direction=f"{ref} higher" if z > 0 else f"{nm} higher"))
    return rows


def wilson(k, n):
    lo, hi = proportion_confint(k, n, alpha=0.05, method="wilson")
    return [float(lo), float(hi)]


def winsorise(x, pct=95):
    cap = float(np.percentile(x, pct))
    return np.minimum(np.asarray(x, float), cap), cap


def location_type(area):
    a = area.lower()
    if "river" in a:
        return "River"
    if "lake" in a or "reservoir" in a:
        return "Lake/reservoir"
    if "port" in a:
        return "Port"
    return "Sea"  # seas, straits, bays and oceans


def accident_group(category):
    """Accident-type group used for the 23 coded Pelra cases and Figure 12."""
    if category in ("Fire", "Explosion"):
        return "Fire/Explosion"
    if "Collision" in category:  # includes 'Capsizing after Collision' (Case 23)
        return "Collision/Contact"
    return category


# ---------------------------------------------------------------- A. accident-level data
raw = pd.read_excel(DATA / "ntsc_maritime_accidents_2007_2022.xlsx", sheet_name="Raw_Data")
raw = raw.rename(columns={"Accident Category": "category", "Analytical Fleet": "fleet", "Final Fatalities": "fat",
                          "Case": "case", "Area": "area", "Year": "year", "Report": "report", "Class": "cas_class"})
raw["year"] = raw["year"].astype(int)
raw["fat"] = raw["fat"].astype(float)
raw["fatal"] = raw["fat"] > 0  # replaces the workbook formula =IF(H>0,"Fatal","Non-fatal")
raw["location"] = raw["area"].map(location_type)
raw["acc_group"] = raw["category"].map(accident_group)
assert len(raw) == 215 and (raw.fleet == "Pelra").sum() == 30
# Three counts corrected against the KNKT reports (dead + missing, crew + passengers); guard against stale copies
corr = pd.read_csv(DATA / "fatality_corrections.csv").set_index("row_id")
assert (raw.set_index("Row ID").loc[corr.index, "fat"] == corr.corrected).all(), "uncorrected raw data"
results["fatality_corrections"] = corr.reset_index().to_dict("records")
acorr = pd.read_csv(DATA / "area_corrections.csv").set_index("row_id")
assert (raw.set_index("Row ID").loc[acorr.index, "area"] == acorr.corrected).all(), "uncorrected area labels"
results["area_corrections"] = acorr.reset_index().to_dict("records")

fleet_n = raw.fleet.value_counts()
retained = [f for f in FLEET_ORDER if fleet_n[f] >= MIN_FLEET_N]
minor = sorted(fleet_n[fleet_n < MIN_FLEET_N].index, key=lambda f: (-fleet_n[f], f))
assert retained == FLEET_ORDER and set(fleet_n.index) == set(FLEET_ORDER) | set(minor)
raw["retained"] = raw.fleet.isin(retained)

fleet_rows = []
for f in FLEET_ORDER + minor:
    d = describe(raw.loc[raw.fleet == f, "fat"])
    fleet_rows.append(dict(fleet=f, retained=f in retained, **d, pct_fatal_ci95=wilson(d["fatal"], d["n"])))
for label, mask in [("Total, seven retained categories", raw.retained), ("Total, all 215 occurrences", raw.fat >= 0)]:
    d = describe(raw.loc[mask, "fat"])
    fleet_rows.append(dict(fleet=label, retained=None, **d, pct_fatal_ci95=wilson(d["fatal"], d["n"])))
fleet_tab = pd.DataFrame(fleet_rows)
fleet_tab.to_csv(OUT / "table3_table10_fleet.csv", index=False)
results["fleet_table"] = fleet_rows
results["minor_categories"] = {f: int(fleet_n[f]) for f in minor}

# Fleet inference on the seven retained categories (n = 208)
groups = [raw.loc[raw.fleet == f, "fat"].to_numpy() for f in FLEET_ORDER]
results["fleet_kruskal"] = kruskal(groups, FLEET_ORDER)
dunn_rows = dunn(groups, FLEET_ORDER, "Pelra", n_contrasts=len(FLEET_ORDER) - 1)
pd.DataFrame(dunn_rows).to_csv(OUT / "dunn_pelra_vs_comparators.csv", index=False)
results["dunn"] = dunn_rows
pelra = raw.loc[raw.fleet == "Pelra", "fat"].to_numpy()
comparators = raw.loc[raw.retained & (raw.fleet != "Pelra"), "fat"].to_numpy()
results["pelra_vs_comparators"] = mwu(pelra, comparators)  # primary: six retained comparator categories
results["pelra_vs_all_nonpelra"] = mwu(pelra, raw.loc[raw.fleet != "Pelra", "fat"].to_numpy())  # sensitivity
passenger = raw.loc[raw.fleet == "Passenger", "fat"].to_numpy()
results["pelra_vs_passenger"] = dict(
    mwu=mwu(pelra, passenger),
    fisher_fatal=dict(zip(["odds_ratio", "p"], map(float, stats.fisher_exact(
        [[int((pelra > 0).sum()), int((pelra == 0).sum())], [int((passenger > 0).sum()), int((passenger == 0).sum())]])))),
    mwu_fatal_only=mwu(pelra[pelra > 0], passenger[passenger > 0]))

# Normality (justification of non-parametric tests)
w_raw = stats.shapiro(pelra)
pelra_w, cap30 = winsorise(pelra)
w_win = stats.shapiro(pelra_w)
results["shapiro_pelra"] = dict(raw=dict(W=float(w_raw.statistic), p=float(w_raw.pvalue), n=len(pelra)),
                                winsorised=dict(W=float(w_win.statistic), p=float(w_win.pvalue), cap=cap30))

# Location types (Figure 12) derived from the Area column
loc = (raw.groupby(["location"]).size()).to_dict()
loc_by = raw.groupby(["acc_group", "fatal", "location"]).size().unstack(fill_value=0)
loc_by.to_csv(OUT / "figure12_location_by_type.csv")
results["location_types"] = {k: int(v) for k, v in loc.items()}
raw[["Row ID", "year", "case", "category", "acc_group", "area", "location", "cas_class", "fleet", "retained", "fat",
     "fatal", "Consequence", "Factor", "report"]].to_csv(OUT / "S1_accident_occurrences.csv", index=False)

# ---------------------------------------------------------------- B. the 30 Pelra cases
key = pd.read_csv(DATA / "coded_case_key.csv")
pel = raw[raw.fleet == "Pelra"].copy()
pel["case_no"] = pel["case"].str.extract(r"(\d+)", expand=False).astype(int)
pel = pel.merge(key, left_on="case_no", right_on="table8_case", how="left")
pel["included"] = pel.coding_case.notna()
assert pel.included.sum() == 23 and (pel.loc[~pel.included, "report"] == "In progress").all()
assert (pel.loc[pel.included, "report"] == "Final").all()
pel["exclusion_reason"] = np.where(pel.included, "", "NTSC report status 'In progress' (not final)")
pel = pel.merge(pd.read_csv(DATA / "pelra_case_attributes.csv"), on="case_no", how="left")
assert pel.classification_group.notna().all()
pel.to_csv(OUT / "S2_pelra_30_cases.csv", index=False,
           columns=["case_no", "year", "category", "acc_group", "area", "location", "cas_class", "fat", "report",
                    "included", "exclusion_reason", "coding_case", "vessel_names", "class_column_J",
                    "classification_group", "knkt_report", "non_wooden_hull", "report_status_cleaned_data"])
inc, exc = pel[pel.included], pel[~pel.included]
results["included_vs_excluded"] = dict(
    mwu=mwu(inc.fat, exc.fat, exact=True),
    excluded_cases=exc.case_no.tolist(), excluded_types=exc.acc_group.value_counts().to_dict(),
    excluded_fatalities=exc.fat.tolist(), included_types=inc.acc_group.value_counts().to_dict())

# Vessel classification within the 30 Pelra cases (groups from column J); empty groups are dropped
CLASS_ORDER = ["Non-Classed Only", "Mixed Class", "Classed Only"]
cls_n = pel.classification_group.value_counts().reindex(CLASS_ORDER, fill_value=0)
names = [nm for nm in CLASS_ORDER if cls_n[nm] > 0]
g = [pel.loc[pel.classification_group == nm, "fat"].to_numpy() for nm in names]
results["classification_groups"] = cls_n.astype(int).to_dict()
results["classification_cases"] = {nm: pel.loc[pel.classification_group == nm, "case_no"].tolist() for nm in names}
if len(names) == 2:
    results["classification_test"] = dict(test="Mann-Whitney", groups=names, **mwu(g[0], g[1], exact=True))
else:
    results["classification_test"] = dict(test="Kruskal-Wallis", **kruskal(g, names))
    results["classification_dunn"] = [r for nm in names for r in dunn(g, names, nm, 3)]

# ---------------------------------------------------------------- C. SHIELD coding (23 cases)
inst = pd.read_csv(DATA / "shield_coding_instances.csv")
assert len(inst) == 114 and inst.table8_case.nunique() == 23
codes = (inst.groupby("code").agg(layer=("layer", "first"), category=("category", "first"), label=("label", "first"))
         .reset_index())
codes["layer_order"] = codes.layer.map(LAYERS.index)
cases23 = inc.sort_values("case_no").case_no.tolist()

presence = pd.crosstab(inst.table8_case, inst.code).reindex(index=cases23, fill_value=0).clip(upper=1)
counts = pd.crosstab(inst.table8_case, inst.code).reindex(index=cases23, fill_value=0)
weights = inst.pivot_table(index="table8_case", columns="code", values="impact_rating", aggfunc="sum",
                           fill_value=0).reindex(index=cases23, fill_value=0)
code_cols = codes.sort_values(["layer_order", "code"]).code.tolist()
for name, m in [("binary", presence), ("instances", counts), ("weights", weights)]:
    m[code_cols].rename_axis("table8_case").to_csv(OUT / f"S3_case_by_code_{name}.csv")

code_tab = codes.set_index("code")
code_tab["cases"] = presence.sum()
code_tab["pct_cases"] = 100 * code_tab.cases / 23
code_tab["instances"] = counts.sum()
code_tab["weight"] = weights.sum()
code_tab["mean_rating"] = code_tab.weight / code_tab.instances
code_tab["median_rating"] = inst.groupby("code").impact_rating.median()
for r in (1, 2, 3, 4):
    code_tab[f"n_rating_{r}"] = inst[inst.impact_rating == r].groupby("code").size().reindex(code_tab.index).fillna(0).astype(int)
code_tab = code_tab.sort_values(["cases", "instances", "weight", "code"], ascending=[False, False, False, True])
code_tab.drop(columns="layer_order").to_csv(OUT / "table12_codes.csv")
assert code_tab.cases.sum() == 109 and code_tab.instances.sum() == 114 and code_tab.weight.sum() == 320

case_layer = inst.groupby(["layer", "table8_case"]).size().unstack(fill_value=0).clip(upper=1)
layer_rows = []
for lay in LAYERS:
    sub = inst[inst.layer == lay]
    layer_rows.append(dict(layer=lay, instances=len(sub), pct_instances=100 * len(sub) / len(inst),
                           cases=int(case_layer.loc[lay].sum()), pct_cases=100 * case_layer.loc[lay].sum() / 23,
                           distinct_codes=sub.code.nunique(), weight=int(sub.impact_rating.sum()),
                           **{f"weight_rating_{r}": int(sub.loc[sub.impact_rating == r, "impact_rating"].sum())
                              for r in (1, 2, 3, 4)},
                           **{f"n_rating_{r}": int((sub.impact_rating == r).sum()) for r in (1, 2, 3, 4)}))
pd.DataFrame(layer_rows).to_csv(OUT / "table1_layers.csv", index=False)
results["layers"] = layer_rows

cat_rows = []
for (lay, cat), sub in inst.groupby(["layer", "category"]):
    cat_rows.append(dict(layer=lay, category=CATEGORY_NAME.get(cat, cat), abbrev=CATEGORY_ABBREV[cat],
                         instances=len(sub), cases=sub.table8_case.nunique(),
                         pct_cases=100 * sub.table8_case.nunique() / 23, weight=int(sub.impact_rating.sum())))
cat_tab = pd.DataFrame(cat_rows)
cat_tab["layer_order"] = cat_tab.layer.map(LAYERS.index)
cat_tab = cat_tab.sort_values(["layer_order", "instances"], ascending=[True, False]).drop(columns="layer_order")
cat_tab.to_csv(OUT / "figure9_17_categories.csv", index=False)
results["categories"] = cat_tab.to_dict("records")

# Table 4: cases with each code by accident-type group
grp = inc.set_index("case_no").acc_group
group_order = ["Capsizing", "Collision/Contact", "Fire/Explosion"]
group_n = grp.value_counts().reindex(group_order)
t4 = pd.DataFrame({g: presence.loc[grp[grp == g].index].sum() for g in group_order})
t4["Total"] = presence.sum()
t4 = t4.loc[code_tab.index]
t4.insert(0, "label", code_tab.label)
t4.to_csv(OUT / "table4_codes_by_accident_type.csv")
results["table4_group_n"] = group_n.astype(int).to_dict()

# ---------------------------------------------------------------- D. RQ3: macro-categories vs fatalities
fat23 = inc.set_index("case_no").fat.reindex(cases23)
fat23_w, cap23 = winsorise(fat23.to_numpy())
fat23_w = pd.Series(fat23_w, index=cases23)
assoc = []
for name, cl in MACRO.items():
    has = presence[[c for c in cl if c in presence.columns]].max(axis=1).astype(bool)
    x, y = fat23[has], fat23[~has]
    r = mwu(x, y, exact=True)
    rw = mwu(fat23_w[has], fat23_w[~has], exact=True)
    assoc.append(dict(macro=name, codes=",".join(cl), cases=sorted(has[has].index.tolist()), **{
        k: r[k] for k in ("n1", "n2", "U_x", "U", "z", "r", "p_asymptotic", "p_exact", "hl_shift", "hl_ci95")},
        median_present=r["x"]["median"], iqr_present=[r["x"]["q1"], r["x"]["q3"]],
        median_absent=r["y"]["median"], iqr_absent=[r["y"]["q1"], r["y"]["q3"]],
        winsor_U_x=rw["U_x"], winsor_p_exact=rw["p_exact"]))
for col in ("p_exact", "p_asymptotic"):
    q = multipletests([a[col] for a in assoc], method="fdr_bh")[1]
    for a, qq in zip(assoc, q):
        a["q_BH_" + col.split("_")[1]] = float(qq)
pd.DataFrame(assoc).to_csv(OUT / "rq3_macro_associations.csv", index=False)
results["rq3_associations"] = assoc
results["rq3_winsor_cap"] = cap23
n_codes, n_inst = presence.sum(axis=1), counts.sum(axis=1)
results["spearman_codes_fatalities"] = {
    "distinct_codes": dict(zip(["rho", "p"], map(float, stats.spearmanr(n_codes, fat23)))),
    "instances": dict(zip(["rho", "p"], map(float, stats.spearmanr(n_inst, fat23))))}
results["shapiro_coded23"] = dict(zip(["W", "p"], map(float, stats.shapiro(fat23))))

# ---------------------------------------------------------------- E. RQ4: intervention priority weights
pk = pd.read_csv(DATA / "packages.csv")
N_STAKEHOLDERS = 15


def package_scores(pk, freq, impact, normalise=False):
    out = []
    for _, row in pk.iterrows():
        cl = row.codes.split(";")
        s = sum(freq[c] * impact[c] for c in cl)
        out.append(s / len(cl) if normalise else s)
    return np.array(out, float)


def rank_desc(v):
    return stats.rankdata(-np.asarray(v), method="min").astype(int)


S_j = 1 + pk.n_support / N_STAKEHOLDERS
F_inst, F_case = code_tab.instances, code_tab.cases
I_mean = code_tab.mean_rating  # F_inst * I_mean == cumulative code weight W_i
I_med = np.ceil(code_tab.median_rating)  # code tier = median rating, half-points rounded up
geo = inst.assign(g=2.0 ** (inst.impact_rating - 1)).groupby("code").g.sum()  # ratings 1,2,3,4 -> 1,2,4,8
pk["codes_display"] = pk.codes.str.replace(";", ", ")
pk["F"] = package_scores(pk, F_inst, pd.Series(1.0, index=F_inst.index))
pk["score"] = package_scores(pk, F_inst, I_mean)
pk["S"] = S_j
pk["PW"] = pk.score * pk.S
pk["rank"] = rank_desc(pk.PW)
pk["case_coverage"] = [int(presence[r.split(";")].max(axis=1).sum()) for r in pk.codes]
sens = {
    "primary": pk.PW,
    "score_only_no_S": pk.score,
    "case_level_tiers": package_scores(pk, F_case, I_med) * S_j,
    "frequency_only": pk.F * S_j,
    "geometric_1_2_4_8": package_scores(pk, geo, pd.Series(1.0, index=geo.index)) * S_j,
    "normalised_per_code": package_scores(pk, F_inst, I_mean, normalise=True) * S_j,
    "case_coverage": pk.case_coverage * S_j,
}
sens_tab = pd.DataFrame({k: rank_desc(v) for k, v in sens.items()})
sens_tab.insert(0, "package", pk.package)
for k in list(sens):
    sens_tab[k + "_value"] = np.round(np.asarray(sens[k], float), 1)
sens_tab.to_csv(OUT / "rq4_sensitivity_ranks.csv", index=False)
rank_agree = {k: dict(zip(["kendall_tau", "p"], map(float, stats.kendalltau(sens_tab.primary, sens_tab[k]))))
              for k in sens if k != "primary"}
pk.to_csv(OUT / "table5_packages.csv", index=False)
results["packages"] = pk.drop(columns=["codes"]).to_dict("records")
results["package_sensitivity_rank_agreement"] = rank_agree
# Intervention audit: for each package and mapped code, the coded cases and evidence excerpts
audit = []
for _, row in pk.iterrows():
    for c in row.codes.split(";"):
        sub = inst[inst.code == c]
        audit.append(dict(package=row.package, code=c, label=code_tab.label[c], cases=int(F_case[c]),
                          instances=int(F_inst[c]), weight=int(code_tab.weight[c]),
                          table8_cases=",".join(map(str, sorted(sub.table8_case.unique()))),
                          evidence=" | ".join(f"C{r.table8_case}: {r.evidence}" for r in sub.itertuples())))
pd.DataFrame(audit).to_csv(OUT / "S4_intervention_audit.csv", index=False)

# ---------------------------------------------------------------- F. sensitivity: non-wooden hulls excluded
# KNKT describes Cases 7 (fibreglass), 26 (fibre) and 28 (steel) as non-wooden; repeat the main analyses without them.
non_wooden = pel.loc[pel.non_wooden_hull.fillna("") != "", "case_no"].tolist()
pel_w = pel[~pel.case_no.isin(non_wooden)]
sub_cases = [c for c in cases23 if c not in non_wooden]
p_sub, i_sub = presence.loc[sub_cases], inst[inst.table8_case.isin(sub_cases)]
f_sub = fat23.loc[sub_cases]
sub_codes = pd.DataFrame({"cases": p_sub.sum(), "instances": counts.loc[sub_cases].sum(),
                          "weight": weights.loc[sub_cases].sum()})
sub_codes["mean_rating"] = (sub_codes.weight / sub_codes.instances).fillna(0)
sub_assoc = []
for name, cl in MACRO.items():
    has = p_sub[cl].max(axis=1).astype(bool)
    rr = mwu(f_sub[has], f_sub[~has], exact=True) if 0 < has.sum() < len(has) else None
    sub_assoc.append(dict(macro=name, n1=int(has.sum()), n2=int((~has).sum()),
                          U=rr["U"] if rr else None, p_exact=rr["p_exact"] if rr else None))
q = multipletests([a["p_exact"] for a in sub_assoc if a["p_exact"] is not None], method="fdr_bh")[1]
for a, qq in zip([a for a in sub_assoc if a["p_exact"] is not None], q):
    a["q_BH_exact"] = float(qq)
pw_sub = package_scores(pk, sub_codes.instances, sub_codes.mean_rating) * S_j
comp_sub = mwu(pel_w.fat.to_numpy(), comparators)
results["sensitivity_non_wooden_excluded"] = dict(
    excluded_cases=non_wooden, n_pelra=len(pel_w), n_coded=len(sub_cases),
    pelra=describe(pel_w.fat), pelra_vs_comparators={k: comp_sub[k] for k in ("U", "U_x", "z", "p_asymptotic", "r")},
    fleet_kruskal={k: v for k, v in kruskal([pel_w.fat.to_numpy()] + groups[1:], FLEET_ORDER).items()
                   if k != "groups"},
    layers={lay: dict(instances=int((i_sub.layer == lay).sum()),
                      pct_instances=100 * (i_sub.layer == lay).mean(),
                      cases=int(i_sub[i_sub.layer == lay].table8_case.nunique())) for lay in LAYERS},
    top_codes={c: dict(cases=int(sub_codes.cases[c]), pct=100 * sub_codes.cases[c] / len(sub_cases))
               for c in ["OR3", "AD1", "AP1", "OC1"]},
    rq3=sub_assoc,
    packages=dict(zip(pk.package, map(float, np.round(pw_sub, 1)))),
    package_ranks=dict(zip(pk.package, map(int, rank_desc(pw_sub)))),
    rank_agreement_tau=float(stats.kendalltau(pk["rank"], rank_desc(pw_sub))[0]))

# ---------------------------------------------------------------- G. other checks
t2 = dict(n=[4, 4, 3, 4], years=[18.5, 21.0, 14.3, 24.8])
results["table2_weighted_mean_experience"] = float(np.dot(t2["n"], t2["years"]) / sum(t2["n"]))

# ---------------------------------------------------------------- H. figures
plt.rcParams.update({"font.family": "serif", "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "savefig.dpi": 300, "savefig.bbox": "tight"})
LAYER_COLOUR = {"Acts": "#2a6fbb", "Preconditions": "#2e8b57", "Operational Leadership": "#e0a020",
                "Organisation": "#8e5bb5"}
FLEET_COLOUR = "#2a6fbb"


def save(fig, name):
    fig.savefig(FIG / f"{name}.png")
    plt.close(fig)


# Figure 4: Pelra vs the six retained comparator categories, broken y-axis
r4 = results["pelra_vs_comparators"]
fig, (top, bot) = plt.subplots(2, 1, sharex=True, figsize=(5.2, 4.2), gridspec_kw=dict(height_ratios=[1, 2.6]))
for ax in (top, bot):
    ax.boxplot([pelra, comparators], widths=0.45, whis=1.5, patch_artist=True,
               boxprops=dict(facecolor="#9ec5e8"), medianprops=dict(color="black"),
               flierprops=dict(marker="o", markersize=3.5, markerfacecolor="none"))
    ax.grid(axis="y", ls=":", lw=0.6)
top.set_ylim(60, 345)
bot.set_ylim(-1, 32)
top.spines["bottom"].set_visible(False)
top.tick_params(bottom=False)
for ax, y in ((top, 0), (bot, 1)):
    ax.plot([0, 0.02], [y - 0.02, y + 0.02], transform=ax.transAxes, color="k", clip_on=False, lw=0.8)
bot.set_xticks([1, 2], [f"Pelra\n(n = {len(pelra)})", f"Six comparator fleet categories\n(n = {len(comparators)})"])
fig.supylabel("Recorded fatalities per occurrence")
top.set_title(f"Mann–Whitney U = {r4['U']:.0f}, p < 0.001, r = {r4['r']:.2f}", fontsize=9, style="italic")
save(fig, "figure4_pelra_vs_comparators")

# Figures 6 and 7: fleet bars (two-word names on two lines so the axis labels do not overlap)
FLEET_LABEL = [f.replace(" ", "\n", 1) for f in FLEET_ORDER]
ft = fleet_tab.set_index("fleet").loc[FLEET_ORDER]
fig, axes = plt.subplots(2, 1, figsize=(6, 5.4))
for ax, col, lab, tag in ((axes[0], "n", "Accident occurrences (n)", "A"),
                          (axes[1], "total", "Total recorded fatalities", "B")):
    b = ax.bar(FLEET_ORDER, ft[col], color=FLEET_COLOUR)
    ax.bar_label(b, fmt="%.0f", fontsize=8)
    ax.set_ylabel(lab)
    ax.set_title(tag, loc="left", fontweight="bold")
    ax.set_xticks(range(len(FLEET_ORDER)), FLEET_LABEL, fontsize=8)
fig.tight_layout()
save(fig, "figure6_fleet_occurrences_fatalities")

fig, axes = plt.subplots(2, 1, figsize=(6, 5.4))
b = axes[0].bar(FLEET_ORDER, ft["mean"], color=FLEET_COLOUR)
axes[0].bar_label(b, fmt="%.2f", fontsize=8)
axes[0].set_ylabel("Mean fatalities per occurrence")
ci = np.array(ft.pct_fatal_ci95.tolist()) * 100
b = axes[1].bar(FLEET_ORDER, ft.pct_fatal, color=FLEET_COLOUR,
                yerr=[ft.pct_fatal - ci[:, 0], ci[:, 1] - ft.pct_fatal], capsize=3, ecolor="black", error_kw=dict(lw=0.8))
axes[1].bar_label(b, labels=[f"{v:.1f}%" for v in ft.pct_fatal], fontsize=8, padding=2)
axes[1].set_ylabel("Occurrences with ≥1 fatality (%)")
axes[1].set_ylim(0, 105)
for ax, tag in zip(axes, "AB"):
    ax.set_title(tag, loc="left", fontweight="bold")
    ax.set_xticks(range(len(FLEET_ORDER)), FLEET_LABEL, fontsize=8)
fig.tight_layout()
save(fig, "figure7_fleet_severity")

# Figure 8: log10(fatalities + 1) by fleet
fig, ax = plt.subplots(figsize=(6.2, 3.4))
vals = [np.log10(g + 1) for g in groups]
ax.boxplot(vals, whis=1.5, widths=0.5, patch_artist=True, boxprops=dict(facecolor="#9ec5e8"),
           medianprops=dict(color="black"), flierprops=dict(marker="o", markersize=3, markerfacecolor="none"))
rng = np.random.default_rng(SEED)
for i, v in enumerate(vals, 1):
    ax.scatter(i + rng.uniform(-0.12, 0.12, len(v)), v, s=5, color="black", alpha=0.35, zorder=3)
ax.set_xticks(range(1, 8), [f"{lab}\n(n = {len(g)})" for lab, g in zip(FLEET_LABEL, groups)], fontsize=7.5)
ax.set_ylabel("log10(fatalities + 1)")
ax.grid(axis="y", ls=":", lw=0.6)
save(fig, "figure8_fleet_log_fatalities")

# Figure 9: coded instances by SHIELD category (cases in parentheses)
fig, ax = plt.subplots(figsize=(7, 3.2))
xs = np.arange(len(cat_tab))
b = ax.bar(xs, cat_tab.instances, color=[LAYER_COLOUR[l] for l in cat_tab.layer])
ax.bar_label(b, labels=[f"{i}\n({c})" for i, c in zip(cat_tab.instances, cat_tab.cases)], fontsize=7)
ax.set_xticks(xs, cat_tab.abbrev)
ax.set_ylabel("Coded instances (n = 114)")
ax.set_ylim(0, cat_tab.instances.max() * 1.25)
ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=c) for c in LAYER_COLOUR.values()], labels=LAYERS,
          ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.12), frameon=False, fontsize=8)
save(fig, "figure9_categories")

# Figure 10: cumulative weight by layer and impact rating
fig, ax = plt.subplots(figsize=(6.2, 3.2))
lt = pd.DataFrame(layer_rows).set_index("layer")
width = 0.2
shades = {4: "#08306b", 3: "#2171b5", 2: "#6baed6", 1: "#c6dbef"}
for j, r in enumerate((4, 3, 2, 1)):
    b = ax.bar(np.arange(4) + (j - 1.5) * width, lt[f"weight_rating_{r}"], width, color=shades[r],
               edgecolor="black", lw=0.4, label=f"Impact rating {r}")
    ax.bar_label(b, fontsize=7)
ax.set_xticks(range(4), [l.replace(" L", "\nL") + f"\n(total = {lt.weight[l]})" for l in LAYERS])
ax.set_ylabel("Cumulative weight (Σ ratings)")
ax.legend(frameon=False, fontsize=8, ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.12))
save(fig, "figure10_layer_weights")

# Figure 12: location types (derived from the Area column)
grp_order = ["Capsizing", "Collision/Contact", "Fire/Explosion", "Grounding", "Other"]
loc_order = ["Sea", "River", "Lake/reservoir", "Port"]
loc_colour = dict(zip(loc_order, ["#2a6fbb", "#7f8c9a", "#d9622b", "#f2b701"]))
fig = plt.figure(figsize=(7, 5.4))
ax = fig.add_axes([0.08, 0.5, 0.9, 0.44])
xpos, labels = [], []
x = 0
for fatal in (True, False):
    for g in grp_order:
        row = loc_by.loc[(g, fatal)] if (g, fatal) in loc_by.index else pd.Series(0, index=loc_order)
        for k, l in enumerate(loc_order):
            v = int(row.get(l, 0))
            bb = ax.bar(x + (k - 1.5) * 0.2, v, 0.2, color=loc_colour[l], label=l if (x == 0) else None)
            if v:
                ax.bar_label(bb, fontsize=6)
        xpos.append(x)
        labels.append(g.replace("/", "/\n"))
        x += 1
    x += 0.5
ax.set_xticks(xpos, labels, fontsize=7)
ymax = int(loc_by.to_numpy().max())
ax.set_ylim(0, ymax * 1.35)
ax.axvline(len(grp_order) - 0.25, color="grey", ls="--", lw=0.7)
ax.text(2, ymax * 1.25, "Fatal", ha="center", fontweight="bold")
ax.text(len(grp_order) + 2.5, ymax * 1.25, "Non-fatal", ha="center", fontweight="bold")
ax.set_ylabel("Accident occurrences")
ax.legend(frameon=False, fontsize=7.5, ncol=4, loc="upper center", bbox_to_anchor=(0.5, 0.88))
ax.set_title("(a)", loc="left", fontweight="bold")
ax2 = fig.add_axes([0.25, 0.02, 0.5, 0.38])
tot = [results["location_types"].get(l, 0) for l in loc_order]
wedges, _ = ax2.pie(tot, colors=[loc_colour[l] for l in loc_order], startangle=90, counterclock=False)
ax2.legend(wedges, [f"{l}: {t} ({100 * t / 215:.1f}%)" for l, t in zip(loc_order, tot)], loc="center left",
           bbox_to_anchor=(1.0, 0.5), frameon=False, fontsize=8)
ax2.set_title("(b)  n = 215", loc="left", fontweight="bold")
save(fig, "figure12_location_types")

# Figures 13-16: code prevalence by layer (cases of 23)
for num, lay in zip((13, 14, 15, 16), LAYERS):
    sub = code_tab[code_tab.layer == lay].sort_values(["cases", "code"], ascending=[False, True])
    fig, ax = plt.subplots(figsize=(6, 0.32 * len(sub) + 0.9))
    b = ax.barh([f"{c} – {l}" for c, l in zip(sub.index, sub.label)][::-1], sub.pct_cases[::-1],
                color=LAYER_COLOUR[lay])
    ax.bar_label(b, labels=[f"{p:.1f}% ({k}/23)" for p, k in zip(sub.pct_cases[::-1], sub.cases[::-1])], fontsize=7,
                 padding=2)
    ax.set_xlim(0, 100)
    ax.set_xlabel("Cases with the code (% of 23 analysed Pelra cases)")
    ax.tick_params(axis="y", labelsize=7.5)
    save(fig, f"figure{num}_{lay.lower().replace(' ', '_')}_codes")

# Figure 17: category prevalence (A) and cumulative weight (B)
c17 = cat_tab.sort_values(["pct_cases", "weight"], ascending=False)
fig, axes = plt.subplots(2, 1, figsize=(7, 5.6))
for ax, col, lab, tag in ((axes[0], "pct_cases", "Cases (% of 23)", "A"),
                          (axes[1], "weight", "Cumulative weight (Σ ratings)", "B")):
    b = ax.bar(c17.category, c17[col], color=[LAYER_COLOUR[l] for l in c17.layer])
    ax.bar_label(b, labels=[f"{v:.1f}" if col == "pct_cases" else f"{v:.0f}" for v in c17[col]], fontsize=7)
    ax.set_ylabel(lab)
    ax.set_title(tag, loc="left", fontweight="bold")
    ax.tick_params(axis="x", labelrotation=55, labelsize=7)
    plt.setp(ax.get_xticklabels(), ha="right")
axes[0].set_ylim(0, 100)
axes[0].legend(handles=[plt.Rectangle((0, 0), 1, 1, color=c) for c in LAYER_COLOUR.values()], labels=LAYERS,
               ncol=4, frameon=False, fontsize=7, loc="upper right")
fig.tight_layout()
save(fig, "figure17_category_prevalence_weight")

# Figure 18: code prevalence and cumulative weight
c18 = code_tab.sort_values(["cases", "weight"], ascending=[True, True])
fig, ax = plt.subplots(figsize=(6.8, 8.2))
yl = [f"{c} – {l}" for c, l in zip(c18.index, c18.label)]
b = ax.barh(yl, c18.pct_cases, color=[LAYER_COLOUR[l] for l in c18.layer])
ax.bar_label(b, labels=[f"{p:.1f}%" for p in c18.pct_cases], fontsize=6.5, padding=2)
ax.set_xlim(0, 100)
ax.set_xlabel("Cases with the code (% of 23)")
ax.tick_params(axis="y", labelsize=7)
ax2 = ax.twiny()
ax2.scatter(c18.weight, yl, color="black", s=12, zorder=3, label="Cumulative weight")
ax2.set_xlim(0, 70)
ax2.set_xlabel("Cumulative weight (Σ ratings, points)")
ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=c) for c in LAYER_COLOUR.values()] +
          [plt.Line2D([], [], marker="o", color="black", ls="")], labels=LAYERS + ["Cumulative weight"],
          frameon=False, fontsize=7, loc="lower right")
save(fig, "figure18_code_prevalence_weight")

# Figure 5: Pelra fatalities by vessel-classification group (broken y-axis for Case 15)
cls_vals = [pel.loc[pel.classification_group == nm, "fat"].to_numpy() for nm in names]
fig, (top, bot) = plt.subplots(2, 1, sharex=True, figsize=(5.2, 4.0), gridspec_kw=dict(height_ratios=[1, 2.6]))
for ax in (top, bot):
    ax.boxplot(cls_vals, widths=0.45, whis=1.5, patch_artist=True, showfliers=False,
               boxprops=dict(facecolor="#9ec5e8"), medianprops=dict(color="black"))
    for i, v in enumerate(cls_vals, 1):
        ax.scatter(i + rng.uniform(-0.1, 0.1, len(v)), v, s=12, color="black", zorder=3)
    ax.grid(axis="y", ls=":", lw=0.6)
top.set_ylim(150, 175)
bot.set_ylim(-1, 32)
top.spines["bottom"].set_visible(False)
top.tick_params(bottom=False)
bot.set_xticks(range(1, len(names) + 1), [f"{nm}\n(n = {len(v)})" for nm, v in zip(names, cls_vals)])
fig.supylabel("Recorded fatalities per occurrence")
ct = results["classification_test"]
top.set_title(f"Mann–Whitney U = {ct['U']:.0f}, exact p = {ct['p_exact']:.2f}" if ct["test"] == "Mann-Whitney"
              else f"Kruskal–Wallis H = {ct['H']:.2f}, p = {ct['p']:.3f}", fontsize=9, style="italic")
save(fig, "figure5_classification")


def box(ax, xy, w, h, title, lines, fc="white", fs=7.5, ls="-"):
    """Box with a bold title and left-aligned lines wrapped to the box width (canvas: 0.75 in per x unit, about 10.5 characters per unit at 7.5 pt)."""
    ax.add_patch(plt.Rectangle(xy, w, h, fc=fc, ec="black", lw=0.8, ls=ls))
    ax.text(xy[0] + w / 2, xy[1] + h - 0.15, title, ha="center", va="top", fontweight="bold", fontsize=fs + 0.5)
    width = int(w * 10.5 * 7.5 / fs) - 1
    wrapped = [w_ for ln in lines for w_ in (textwrap.wrap(ln, width, subsequent_indent="   ") or [""])]
    ax.text(xy[0] + 0.15, xy[1] + h - 0.5, "\n".join(wrapped), ha="left", va="top", fontsize=fs, linespacing=1.4)


def arrow(ax, a, b, ls="-"):
    ax.annotate("", xy=b, xytext=a, arrowprops=dict(arrowstyle="-|>", lw=0.9, ls=ls, color="black"))


def canvas(height):
    fig, ax = plt.subplots(figsize=(9, height * 0.75))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, height)
    ax.axis("off")
    return fig, ax


# Figure 1: study design; the stakeholder pathway feeds the priority weights
fig, ax = canvas(9.6)
box(ax, (0.1, 6.7), 3.8, 2.8, "Documentary analysis (primary)",
    [f"NTSC (KNKT) investigations, 2007–2022: n = {len(raw)} accident occurrences",
     f"Fleet comparison: {len(retained)} commercial categories with ≥{MIN_FLEET_N} occurrences "
     f"(n = {int(raw.retained.sum())})",
     f"Pelra subset: n = {len(pel)}", f"SHIELD coding of final reports: n = {len(inc)}"])
box(ax, (4.1, 6.7), 3.8, 2.8, "Literature and standards",
    ["SHIELD taxonomy, four layers [20]", "HFACS-based maritime accident studies",
     "IMO casualty severity classes (VSMC, SMC)"])
box(ax, (8.1, 6.7), 3.8, 2.8, "Stakeholder consultation",
    [f"n = {N_STAKEHOLDERS} purposively selected experts: NTSC investigators (4), MoT/Syahbandar officials (4), "
     "BKI surveyors (3), Pelra operators and shipowners (4)",
     "Support votes for the 10 intervention packages (n_support)"])
box(ax, (0.3, 2.9), 7.4, 3.2, "Analysis",
    ["RQ1  Fleet severity: Kruskal–Wallis, Dunn (Bonferroni), Mann–Whitney U",
     "RQ2  SHIELD codes and layers: 23 × 32 incidence matrix, impact ratings 1–4",
     "RQ3  Factor–fatality associations: exact Mann–Whitney U, Benjamini–Hochberg",
     "RQ4  Priority weights PWj = Sj × Σ Fi·Ii (Equation (1))"])
box(ax, (8.3, 3.4), 3.4, 2.0, "Stakeholder modifier",
    ["Sj = 1 + n_support / 15 (Equation (2))", "Multiplies the accident-derived score in RQ4"], fc="#eef4fb")
box(ax, (1.8, 0.2), 8.4, 1.9, "Results and conclusions",
    ["Casualty severity conditional on an investigated occurrence; recurrent SHIELD codes and layers; "
     "factor–fatality associations; prioritised intervention packages and NTSC recommendation cross-walk"])
arrow(ax, (2.0, 6.7), (2.6, 6.1))
arrow(ax, (6.0, 6.7), (5.4, 6.1))
arrow(ax, (10.0, 6.7), (10.0, 5.4))
arrow(ax, (8.3, 4.1), (7.7, 3.4))
arrow(ax, (4.0, 2.9), (4.8, 2.1))
save(fig, "figure1_study_design")

# Figure 2: case-selection flow
fig, ax = canvas(9.0)
minor_txt = ", ".join(f"{f.lower()} {n}" for f, n in results["minor_categories"].items())
exc_txt = ", ".join(map(str, exc.case_no))
gn = results["table4_group_n"]
box(ax, (3.3, 7.6), 5.4, 1.2, "NTSC maritime accident occurrences, 2007–2022", [f"n = {len(raw)}"])
box(ax, (0.1, 5.0), 5.4, 1.8, "Fleet comparison (RQ1)",
    [f"{len(retained)} commercial fleet categories with ≥{MIN_FLEET_N} occurrences: n = {int(raw.retained.sum())} "
     f"(Pelra n = {len(pel)})"])
box(ax, (0.1, 2.6), 5.4, 1.8, "Descriptive inventory only",
    [f"{len(minor)} categories with <{MIN_FLEET_N} occurrences, n = {int((~raw.retained).sum())}: {minor_txt}"],
    fc="#f3f3f3")
box(ax, (6.5, 5.0), 5.4, 1.8, "Pelra vessel subset", [f"n = {len(pel)}"])
box(ax, (8.7, 2.6), 3.2, 1.8, "Excluded from SHIELD coding",
    [f"n = {len(exc)}: NTSC report status 'In progress' (not final)", f"Cases {exc_txt}"], fc="#f3f3f3", fs=7)
box(ax, (6.5, 0.2), 5.4, 1.8, "SHIELD-coded final reports (RQ2–RQ4)",
    [f"n = {len(inc)}: capsizing {gn['Capsizing']}, collision/contact {gn['Collision/Contact']}, "
     f"fire/explosion {gn['Fire/Explosion']}"])
arrow(ax, (5.0, 7.6), (2.8, 6.8))
arrow(ax, (7.0, 7.6), (9.2, 6.8))
arrow(ax, (2.8, 5.0), (2.8, 4.4), ls="--")
arrow(ax, (7.4, 5.0), (7.4, 2.0))
arrow(ax, (7.4, 3.5), (8.7, 3.5), ls="--")
save(fig, "figure2_case_selection")


# ---------------------------------------------------------------- I. save
def clean(o):
    if isinstance(o, dict):
        return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating, float)):
        return None if math.isnan(o) else float(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


(OUT / "results.json").write_text(json.dumps(clean(results), indent=1, ensure_ascii=False))
print("done:", len(list(OUT.iterdir())), "outputs,", len(list(FIG.iterdir())), "figures")
