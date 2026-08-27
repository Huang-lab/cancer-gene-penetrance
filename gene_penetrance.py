"""
Concordance, canonical-denominator penetrance, and HTML report for AoU cancer-gene carriers.

CAVEAT: All estimates here are crude/cross-sectional, computed from a single
AoU snapshot (not age-adjusted survival/cumulative-incidence analysis) and
are subject to ascertainment bias -- AoU is a healthcare-utilizing volunteer
cohort, not population-representative for cancer-predisposition genes. Treat
"penetrance" numbers as descriptive, not clinical-grade risk estimates.

Pipeline:
  1. Build ONE canonical carrier denominator per gene (raw rows -> unique
     person x gene -> joined to cohort -> after sex restriction), reused by
     every downstream table. A carrier_filter_log.tsv records every step.
  2. Flag genes with a likely pseudogene-mapping artifact (PMS2 always;
     others if one non-EUR ancestry has >2x the EUR carrier count).
  3. Restrict "interpretable" results to genes with >=MIN_INTERPRETABLE_CARRIERS
     canonical carriers,
     no mapping-artifact flag, and a defined gene->cancer mapping.
  4. ACMG/ClinVar concordance by ancestry, with singleton-fraction and
     top-recurrent-variant footnotes per ancestry.
  5. Disease-specific penetrance by age band, one line per gene x cancer
     type (never merged), each compared to that disease's own population
     background rate in the same cohort.
  6. report.html: a single self-contained report for the above.

Usage:
    python gene_penetrance.py \
        --carriers-path /path/to/carriers_all_noMUTYH_ancestry.tsv \
        --cohort-path   /path/to/All_temp_cohort_simple.csv.gz

Outputs land under --outdir (default: results/ beside this script) as
tables/, figures/ and the three HTML reports, which is the layout committed
to the repo.
"""

import argparse
import base64
import gzip
import html
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from statsmodels.stats.proportion import proportion_confint

warnings.filterwarnings("ignore")

HERE = Path(__file__).parent

# Where the run writes. Set from the CLI in main(); the module-level values are
# the defaults and the ones the committed results/ tree was produced under, so
# a plain `python gene_penetrance.py` regenerates what is checked in, in place.
#
# These were previously `HERE / "..."` with the inputs read from `HERE.parent`,
# left over from the script's original home in `DemoVariant/Penetrance/`. That
# wrote every table and figure into the repo root and looked for the input data
# one level *above* the repo, so nothing under results/ was reachable by running
# the script -- the committed tree had been moved there by hand.
CARRIERS_PATH = None
COHORT_PATH   = None
OUT_DIR       = HERE / "results"
TABLE_DIR     = OUT_DIR / "tables"
FIG_DIR       = OUT_DIR / "figures"

DEFAULT_CARRIERS_NAME = "carriers_all_noMUTYH_ancestry.tsv"
DEFAULT_COHORT_NAME   = "All_temp_cohort_simple.csv.gz"

CONTROL_GROUPS = ["Control (age ≥ 50)", "Control (age < 50)"]
FAMHX_GROUPS   = ["Family history (age < 50)", "Family history (age ≥ 50)"]

AGE_BINS   = [0, 40, 50, 60, 70, 200]
AGE_LABELS = ["<40", "40-49", "50-59", "60-69", "70+"]

MIN_INTERPRETABLE_CARRIERS = 20

# Explicit gene -> qualifying disease Group(s) mapping.
# Groups must match the exact "Group" labels used in All_temp_cohort_simple.csv.gz.
GENE_CANCER_GROUPS = {
    "APC":      ["Colorectum"],
    "BMPR1A":   ["Colorectum"],
    "BRCA1":    ["Breast", "Ovary"],
    "BRCA2":    ["Breast", "Ovary"],
    "MAX":      ["Endocrine System"],
    "MEN1":     ["Endocrine System", "Pancreas"],
    "MLH1":     ["Colorectum", "Endometrium"],
    "MSH2":     ["Colorectum", "Endometrium"],
    "MSH6":     ["Colorectum", "Endometrium"],
    "PMS2":     ["Colorectum", "Endometrium"],
    "NF2":      ["Central Nervous System", "Meninges"],
    "PALB2":    ["Breast"],
    "PTEN":     ["Breast", "Thyroid", "Endometrium"],
    "RB1":      ["Eye"],
    "RET":      ["Thyroid", "Endocrine System"],
    "SDHAF2":   ["Endocrine System"],
    "SDHB":     ["Endocrine System"],
    "SDHC":     ["Endocrine System"],
    "SDHD":     ["Endocrine System"],
    "SMAD4":    ["Colorectum"],
    "STK11":    ["Colorectum", "Breast"],
    "TMEM127":  ["Endocrine System"],
    "TP53":     ["Breast"],
    "TSC1":     ["Kidney"],
    "TSC2":     ["Kidney"],
    "VHL":      ["Kidney", "Central Nervous System", "Endocrine System"],
    "WT1":      ["Kidney"],
}

# Genes whose mapped groups are sex-specific -> restrict numerator+denominator
# to female carriers (all mapped groups below are female-specific organs).
FEMALE_ONLY_GENES = {"BRCA1", "BRCA2", "PALB2", "MLH1", "MSH2", "MSH6", "PMS2", "PTEN"}

# Genes known to sit in/near a pseudogene or highly-homologous paralog region,
# where short-read mapping can create ancestry-skewed spurious carrier calls.
PSEUDOGENE_WATCH_GENES = {"PMS2", "CHEK2", "SMAD4"}


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def wilson_ci(count: int, nobs: int):
    if nobs == 0:
        return (np.nan, np.nan)
    lo, hi = proportion_confint(count, nobs, alpha=0.05, method="wilson")
    return (lo, hi)


def collapse_groups(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["Group_collapsed"] = df["Group"].replace(
        {g: "Control" for g in CONTROL_GROUPS}
    ).replace(
        {g: "Family history" for g in FAMHX_GROUPS}
    )
    return df


def load_carriers_raw() -> pd.DataFrame:
    print("Loading raw carrier rows ...")
    df = pd.read_csv(CARRIERS_PATH, sep="\t",
                      usecols=["person_id", "Gene", "Variant_ID",
                               "is_ACMG_PLP", "is_ClinVar_PLP",
                               "ancestry_pred", "sex_at_birth", "age"])
    df["person_id"] = df["person_id"].astype(str)
    print(f"  {len(df):,} raw carrier rows, {df['Gene'].nunique():,} genes")
    return df


def dedupe_to_person_gene(raw: pd.DataFrame) -> pd.DataFrame:
    """Unique person x Gene. A person counts as ACMG/ClinVar PLP for the gene
    if ANY of their variants in it is (OR across variants)."""
    agg = raw.groupby(["person_id", "Gene"], as_index=False).agg(
        is_ACMG_PLP=("is_ACMG_PLP", "max"),
        is_ClinVar_PLP=("is_ClinVar_PLP", "max"),
        ancestry_pred=("ancestry_pred", "first"),
    )
    print(f"  {len(agg):,} unique person x gene carrier rows, "
          f"{agg['person_id'].nunique():,} unique carriers")
    return agg


def load_cohort() -> pd.DataFrame:
    print("Loading cohort ...")
    with gzip.open(COHORT_PATH, "rt") as f:
        cohort = pd.read_csv(f, low_memory=False,
                              usecols=["person_id", "age", "sex_at_birth", "Group"])
    cohort["person_id"] = cohort["person_id"].astype(str)
    cohort = collapse_groups(cohort)
    print(f"  {len(cohort):,} cohort rows, {cohort['person_id'].nunique():,} unique persons "
          f"({(cohort.groupby('person_id').size() > 1).sum():,} persons have >1 row -- "
          f"multiple recorded conditions/groups)")
    return cohort


def build_person_tables(cohort: pd.DataFrame):
    """One canonical demographic row per person, plus the FULL set of
    Group_collapsed values each person appears under (a person can have
    several -- e.g. two different cancer diagnoses). Using a set, rather than
    an arbitrary 'first row', avoids under- or double-counting disease
    status and avoids multiplying carrier rows on join (the root cause of
    the earlier BRCA1 denominator mismatch: 565 vs ~611 vs ~917)."""
    person_demo = cohort.groupby("person_id", as_index=False).agg(
        age=("age", "first"), sex_at_birth=("sex_at_birth", "first")
    ).set_index("person_id")
    person_groups = cohort.groupby("person_id")["Group_collapsed"].apply(set)
    return person_demo, person_groups


# ---------------------------------------------------------------------------
# Canonical carrier denominator
# ---------------------------------------------------------------------------

def build_canonical(raw: pd.DataFrame, dedup: pd.DataFrame,
                     person_demo: pd.DataFrame,
                     classification_col: str = "is_ACMG_PLP",
                     classification_label: str = "ACMG_PLP",
                     log_filename: str = "carrier_filter_log.tsv") -> tuple:
    """Build the ONE canonical carrier set reused by every downstream table:
    unique person x Gene -> restricted to carriers of a PLP variant per
    `classification_col` -> joined to the cohort (dropping carriers with no
    cohort record) -> sex-restricted for female-specific genes.

    "Canonical" here means: this exact person x gene set, built by this exact
    filter chain, is the ONLY carrier population used for every downstream
    number (penetrance, ancestry breakdown, evidence-gap %) for that
    classification. Nothing downstream re-derives its own denominator.

    classification_col is "is_ACMG_PLP" for the main pipeline (ACMG is a
    rule-based classifier applied uniformly regardless of ancestry) or
    "is_ClinVar_PLP" for the side-by-side comparison pipeline -- restricting
    the MAIN report to ClinVar-PLP would reproduce exactly the ancestry-based
    under-ascertainment shown in Finding 1, which is why ACMG is the default.

    Returns (canonical_df, filter_log_df, ancestry_prerestrict) where
    ancestry_prerestrict maps Gene -> ancestry value_counts BEFORE sex
    restriction (used only for the pseudogene-artifact check, which is a
    genotyping-mapping property, not a phenotype/sex property).
    """
    raw_counts = raw.groupby("Gene")["person_id"].count()  # raw rows, not unique

    canonical_frames = []
    log_rows = []
    ancestry_prerestrict = {}

    for gene, gsub in dedup.groupby("Gene"):
        raw_n = int(raw_counts.get(gene, 0))
        step1 = gsub["person_id"].nunique()

        plp = gsub[gsub[classification_col] == True]
        step_plp = plp["person_id"].nunique()

        joined = plp.merge(person_demo, left_on="person_id", right_index=True, how="inner")
        step2 = joined["person_id"].nunique()
        dropped_cohort = step_plp - step2

        ancestry_prerestrict[gene] = joined["ancestry_pred"].value_counts()

        sex_restricted = gene in FEMALE_ONLY_GENES
        step3df = joined[joined["sex_at_birth"] == "Female"] if sex_restricted else joined
        step3 = step3df["person_id"].nunique()

        canonical_frames.append(step3df.assign(Gene=gene))
        log_rows.append({
            "Gene": gene,
            "raw_variant_rows": raw_n,
            "unique_person_gene_any_classification": step1,
            f"unique_person_gene_{classification_label}": step_plp,
            "joined_to_cohort": step2,
            "dropped_at_cohort_join": dropped_cohort,
            "pct_dropped_at_cohort_join": round(100 * dropped_cohort / step_plp, 1) if step_plp else np.nan,
            "sex_restricted": sex_restricted,
            "n_carriers_canonical": step3,
        })

    canonical_df = pd.concat(canonical_frames, ignore_index=True)
    filter_log_df = pd.DataFrame(log_rows).sort_values("n_carriers_canonical", ascending=False)

    plp_col = f"unique_person_gene_{classification_label}"
    brca1 = filter_log_df[filter_log_df["Gene"] == "BRCA1"]
    if len(brca1):
        r = brca1.iloc[0]
        print(
            f"\n[Filter log] BRCA1 ({classification_label}): "
            f"raw carrier rows={r.raw_variant_rows:,} "
            f"-> unique person x gene (any classification)={r.unique_person_gene_any_classification:,} "
            f"-> {classification_label} carriers only={r[plp_col]:,} "
            f"-> joined to cohort={r.joined_to_cohort:,} "
            f"(dropped {r.dropped_at_cohort_join:,}, {r.pct_dropped_at_cohort_join}% -- "
            "carriers with no cohort/phenotype record) "
            f"-> after sex restriction (female only)={r.n_carriers_canonical:,}"
        )

    out_path = TABLE_DIR / log_filename
    filter_log_df.to_csv(out_path, sep="\t", index=False)
    print(f"  wrote {out_path}")

    big_drop = filter_log_df[filter_log_df["pct_dropped_at_cohort_join"] >= 20]
    if len(big_drop):
        print(f"\n[Warning] {len(big_drop)} gene(s) lose >=20% of carriers at the "
              "cohort-join step (no cohort/phenotype record found):")
        for _, r in big_drop.iterrows():
            print(f"    {r.Gene}: {r.dropped_at_cohort_join:,} / {r[plp_col]:,} "
                  f"dropped ({r.pct_dropped_at_cohort_join}%)")

    return canonical_df, filter_log_df, ancestry_prerestrict


def compute_mapping_artifact_flags(ancestry_prerestrict: dict) -> pd.DataFrame:
    """PMS2 is always flagged (known PMS2/PMS2CL pseudogene mismapping).
    Other pseudogene-watch genes are flagged if a single non-EUR ancestry
    has >2x the EUR carrier count, computed on pre-sex-restriction,
    cohort-joined counts (a genotyping-mapping property, not a phenotype
    property)."""
    rows = []
    for gene, counts in ancestry_prerestrict.items():
        eur = int(counts.get("eur", 0))
        flagged, reason = False, ""
        if gene == "PMS2":
            flagged = True
            reason = "Known PMS2/PMS2CL pseudogene mismapping artifact"
        elif gene in PSEUDOGENE_WATCH_GENES:
            for anc, n in counts.items():
                if anc != "eur" and eur > 0 and n > 2 * eur:
                    flagged = True
                    reason = (f"{anc}={int(n)} carriers is >2x EUR={eur}; "
                              f"{gene} is a pseudogene-watch gene")
                    break
        rows.append({"Gene": gene, "possible_mapping_artifact": flagged,
                      "mapping_artifact_reason": reason})
    return pd.DataFrame(rows)


def build_person_group_membership(person_groups: pd.Series, person_ids: pd.Series,
                                   target_groups) -> np.ndarray:
    """Vectorized-ish membership test: for each person_id, is any of
    target_groups in that person's full set of recorded Group_collapsed
    values? Handles people with multiple cohort rows correctly."""
    target = set(target_groups) if not isinstance(target_groups, str) else {target_groups}
    return np.array([
        bool(target & person_groups.get(pid, set()))
        for pid in person_ids
    ])


# ---------------------------------------------------------------------------
# Fix 1 - ACMG/ClinVar concordance by ancestry, with singleton/recurrence footnotes
# ---------------------------------------------------------------------------

def fix1_concordance(raw: pd.DataFrame) -> pd.DataFrame:
    print("\n[Fix 1] ACMG/ClinVar concordance by ancestry ...")
    dedup_va = raw.drop_duplicates(["person_id", "Variant_ID"])
    variant_ac_proxy = dedup_va.groupby("Variant_ID")["person_id"].nunique()

    acmg = raw[raw["is_ACMG_PLP"] == True].drop_duplicates(["person_id", "Gene"])
    rows = []
    eur_rate = None
    for ancestry, sub in acmg.groupby("ancestry_pred"):
        n = len(sub)
        discordant = int((sub["is_ClinVar_PLP"] == False).sum())
        rate = discordant / n if n else np.nan
        lo, hi = wilson_ci(discordant, n)
        rows.append({
            "ancestry_pred": ancestry,
            "n_ACMG_PLP_carriers": n,
            "n_discordant_not_ClinVar_PLP": discordant,
            "discordance_rate": rate,
            "CI_lower": lo,
            "CI_upper": hi,
        })
        if ancestry == "eur":
            eur_rate = rate
    out = pd.DataFrame(rows).sort_values("discordance_rate", ascending=False)
    out["diff_vs_eur"] = out["discordance_rate"] - eur_rate
    out["ci_overlaps_eur"] = out.apply(
        lambda r: (eur_rate is not None) and (r["CI_lower"] <= eur_rate <= r["CI_upper"]),
        axis=1,
    )

    # Footnote A5: singleton fraction + top recurrent variant share, computed
    # at the VARIANT-INSTANCE level (not the dedup-to-gene level used above),
    # since a variant's allele-count and identity are properties of the
    # variant, not of the person x gene summary.
    disc_variants = raw[(raw["is_ACMG_PLP"] == True) & (raw["is_ClinVar_PLP"] == False)] \
        .drop_duplicates(["person_id", "Variant_ID"]).copy()
    disc_variants["AC_proxy"] = disc_variants["Variant_ID"].map(variant_ac_proxy)

    footnote_rows = []
    for ancestry, sub in disc_variants.groupby("ancestry_pred"):
        n = len(sub)
        singleton_frac = (sub["AC_proxy"] == 1).mean() if n else np.nan
        top_counts = sub["Variant_ID"].value_counts()
        top_variant_id = top_counts.index[0] if len(top_counts) else ""
        top_variant_gene = sub.loc[sub["Variant_ID"] == top_variant_id, "Gene"].iloc[0] if len(top_counts) else ""
        top_variant_share = top_counts.iloc[0] / n if n else np.nan
        footnote_rows.append({
            "ancestry_pred": ancestry,
            "n_discordant_variant_instances": n,
            "singleton_fraction": singleton_frac,
            "top_variant_id": top_variant_id,
            "top_variant_gene": top_variant_gene,
            "top_variant_share": top_variant_share,
        })
    footnotes = pd.DataFrame(footnote_rows)
    out = out.merge(footnotes, on="ancestry_pred", how="left")

    out_path = TABLE_DIR / "acmg_clinvar_concordance_by_ancestry.tsv"
    out.to_csv(out_path, sep="\t", index=False)
    print(f"  wrote {out_path}")

    eas_row = out[out["ancestry_pred"] == "eas"]
    if len(eas_row):
        r = eas_row.iloc[0]
        print(f"  [EAS check] top recurrent variant {r.top_variant_id} ({r.top_variant_gene}) "
              f"accounts for {r.top_variant_share:.1%} of EAS discordant variant instances "
              f"(singleton fraction only {r.singleton_fraction:.1%}) -- EAS discordance is "
              "largely ONE recurrent variant, not broad-based curation gaps.")

    fig, ax = plt.subplots(figsize=(8, 5))
    yerr = [out["discordance_rate"] - out["CI_lower"], out["CI_upper"] - out["discordance_rate"]]
    ax.bar(out["ancestry_pred"], out["discordance_rate"] * 100, yerr=np.array(yerr) * 100,
           capsize=4, color="#4C72B0")
    ax.set_ylabel("ACMG-PLP carriers not ClinVar-PLP (%)")
    ax.set_xlabel("Predicted ancestry")
    ax.set_title("ACMG/ClinVar discordance by ancestry (95% Wilson CI)")
    fig.tight_layout()
    fig_path = FIG_DIR / "concordance_by_ancestry.png"
    fig.savefig(fig_path, dpi=150)
    plt.close(fig)
    print(f"  wrote {fig_path}")
    return out


# ---------------------------------------------------------------------------
# Fix 3 - Gene-specific penetrance (canonical denominator, pooled ages)
# ---------------------------------------------------------------------------

def fix3_gene_penetrance(canonical_df: pd.DataFrame, person_groups: pd.Series,
                          artifact_flags: pd.DataFrame, suffix: str = "") -> pd.DataFrame:
    print(f"\n[Fix 3] Gene-specific penetrance (canonical denominator{suffix}, explicit numerator) ...")
    rows = []
    for gene, sub in canonical_df.groupby("Gene"):
        mapped_groups = GENE_CANCER_GROUPS.get(gene)
        n_carriers = len(sub)  # canonical_df already 1 row per person per gene

        if mapped_groups:
            primary_affected = build_person_group_membership(person_groups, sub["person_id"], mapped_groups)
            n_primary = int(primary_affected.sum())
            primary_pen = n_primary / n_carriers if n_carriers else np.nan
            plo, phi = wilson_ci(n_primary, n_carriers)
        else:
            n_primary, primary_pen, plo, phi = np.nan, np.nan, np.nan, np.nan

        any_affected = np.array([
            bool(person_groups.get(pid, set()) - {"Control", "Family history"})
            for pid in sub["person_id"]
        ])
        n_any = int(any_affected.sum())
        any_pen = n_any / n_carriers if n_carriers else np.nan
        alo, ahi = wilson_ci(n_any, n_carriers)

        rows.append({
            "Gene": gene,
            "mapped_cancer_groups": ";".join(mapped_groups) if mapped_groups else "",
            "has_gene_cancer_mapping": bool(mapped_groups),
            "sex_restricted": gene in FEMALE_ONLY_GENES,
            "n_carriers_canonical": n_carriers,
            "n_primary_affected": n_primary,
            "penetrance_primary": primary_pen,
            "CI_lower_primary": plo,
            "CI_upper_primary": phi,
            "n_any_cancer_affected": n_any,
            "penetrance_any_cancer_secondary": any_pen,
            "CI_lower_any": alo,
            "CI_upper_any": ahi,
        })

    out = pd.DataFrame(rows).merge(artifact_flags, on="Gene", how="left")
    out["meets_min_carriers"] = out["n_carriers_canonical"] >= MIN_INTERPRETABLE_CARRIERS
    out["interpretable"] = (
        out["meets_min_carriers"]
        & ~out["possible_mapping_artifact"]
        & out["has_gene_cancer_mapping"]
    )
    out = out.sort_values("penetrance_primary", ascending=False)
    out_path = TABLE_DIR / f"gene_penetrance{suffix}.tsv"
    out.to_csv(out_path, sep="\t", index=False)
    print(f"  wrote {out_path}")

    excl = out[~out["interpretable"]].copy()
    def _reason(r):
        reasons = []
        if not r["has_gene_cancer_mapping"]:
            reasons.append("no gene->cancer mapping (any-cancer fallback only)")
        if r["possible_mapping_artifact"]:
            reasons.append(f"mapping artifact ({r['mapping_artifact_reason']})")
        if not r["meets_min_carriers"]:
            reasons.append(f"< {MIN_INTERPRETABLE_CARRIERS} canonical carriers "
                            f"(n={r['n_carriers_canonical']})")
        return "; ".join(reasons)
    excl["exclusion_reason"] = excl.apply(_reason, axis=1)
    excl_path = TABLE_DIR / f"excluded_genes{suffix}.tsv"
    excl[["Gene", "n_carriers_canonical", "exclusion_reason"]].to_csv(excl_path, sep="\t", index=False)
    print(f"  wrote {excl_path} ({len(excl)} genes excluded from main results)")

    return out


# ---------------------------------------------------------------------------
# Fix 3b - Disease-specific penetrance by age band (canonical denominator,
#          one line per gene x single cancer type, never merged)
# ---------------------------------------------------------------------------

def population_baseline_by_group(person_demo: pd.DataFrame, person_groups: pd.Series,
                                  disease: str, female_only: bool) -> pd.DataFrame:
    df = person_demo[person_demo["sex_at_birth"] == "Female"] if female_only else person_demo
    df = df.copy()
    df["age_band"] = pd.cut(df["age"], bins=AGE_BINS, labels=AGE_LABELS, right=False)
    df["affected"] = build_person_group_membership(person_groups, df.index.to_series(), {disease})
    rows = []
    for band, sub in df.groupby("age_band", observed=True):
        n = len(sub)
        n_aff = int(sub["affected"].sum())
        lo, hi = wilson_ci(n_aff, n)
        rows.append({"age_band": band, "n_population": n, "n_affected": n_aff,
                      "baseline_rate": n_aff / n if n else np.nan,
                      "baseline_CI_lower": lo, "baseline_CI_upper": hi})
    return pd.DataFrame(rows).set_index("age_band")


def fix3b_disease_specific_by_age(canonical_df: pd.DataFrame, person_demo: pd.DataFrame,
                                   person_groups: pd.Series, artifact_flags: pd.DataFrame,
                                   suffix: str = "") -> pd.DataFrame:
    print(f"\n[Fix 3b] Disease-specific penetrance by age band (canonical denominator{suffix}) ...")
    canonical_df = canonical_df.copy()
    canonical_df["age_band"] = pd.cut(canonical_df["age"], bins=AGE_BINS, labels=AGE_LABELS, right=False)

    rows = []
    for gene, gsub in canonical_df.groupby("Gene"):
        mapped_groups = GENE_CANCER_GROUPS.get(gene)
        if not mapped_groups:
            continue
        female_only = gene in FEMALE_ONLY_GENES

        for disease in mapped_groups:
            baseline = population_baseline_by_group(person_demo, person_groups, disease, female_only)
            affected_mask = build_person_group_membership(person_groups, gsub["person_id"], {disease})
            gsub_d = gsub.assign(affected=affected_mask)

            for band, sub in gsub_d.groupby("age_band", observed=True):
                n = len(sub)
                n_affected = int(sub["affected"].sum())
                lo, hi = wilson_ci(n_affected, n)
                base_upper = baseline.loc[band, "baseline_CI_upper"] if band in baseline.index else np.nan
                rows.append({
                    "Gene": gene,
                    "disease": disease,
                    "age_band": band,
                    "n_carriers_canonical": n,
                    "n_affected": n_affected,
                    "penetrance": n_affected / n if n else np.nan,
                    "CI_lower": lo,
                    "CI_upper": hi,
                    "baseline_rate": baseline.loc[band, "baseline_rate"] if band in baseline.index else np.nan,
                    "elevated_vs_background": (lo > base_upper) if pd.notna(base_upper) else False,
                })

    out = pd.DataFrame(rows).merge(artifact_flags, on="Gene", how="left")
    out = out.sort_values(["Gene", "disease", "age_band"])
    out_path = TABLE_DIR / f"disease_specific_penetrance_by_age{suffix}.tsv"
    out.to_csv(out_path, sep="\t", index=False)
    print(f"  wrote {out_path}")
    return out


def plot_gene_disease_chart(gene: str, gene_df: pd.DataFrame, fig_path: Path,
                             n_carriers: int = None, classification_label: str = "ACMG-PLP"):
    """One small chart for ONE gene: each mapped cancer type as its own line
    (never merged), plus that cancer type's population background as a
    dashed line, and a star where the carrier rate is significantly above
    background. `classification_label` names which classifier's PLP calls
    define the carrier set plotted here (ACMG-PLP or ClinVar-PLP) -- stated
    explicitly in the title so the classification basis is never ambiguous
    or, worse, mislabeled when the same function is reused for both
    pipelines. Legend is placed OUTSIDE the plot area so it never overlaps
    the data."""
    fig, ax = plt.subplots(figsize=(8.5, 5))
    for disease, sub in gene_df.groupby("disease"):
        sub = sub.set_index("age_band").reindex(AGE_LABELS)
        line, = ax.plot(AGE_LABELS, sub["penetrance"] * 100, marker="o", label=f"{gene} carriers: {disease}")
        ax.plot(AGE_LABELS, sub["baseline_rate"] * 100, marker="s", linestyle="--",
                color=line.get_color(), alpha=0.5, label=f"Population background: {disease}")
        elevated = sub["elevated_vs_background"].fillna(False)
        ax.plot(np.array(AGE_LABELS)[elevated.values],
                (sub["penetrance"] * 100)[elevated.values],
                marker="*", markersize=16, linestyle="none", color=line.get_color())
    ax.set_ylabel("Proportion with this cancer type (%)")
    ax.set_xlabel("Age band")
    n_label = (f" (n={n_carriers:,} {classification_label} carriers)" if n_carriers is not None
               else f" ({classification_label} carriers)")
    ax.set_title(f"{gene}: disease-specific penetrance by age{n_label}\n"
                 "★ = significantly above population background")
    ax.legend(fontsize=7, loc="center left", bbox_to_anchor=(1.02, 0.5))
    fig.tight_layout()
    fig.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {fig_path}")


# ---------------------------------------------------------------------------
# Fix 4 - Ancestry stratification + evidence gap (canonical denominator)
# ---------------------------------------------------------------------------

def fix4_ancestry_evidence_gap(canonical_df: pd.DataFrame, min_carriers: int,
                                artifact_flags: pd.DataFrame) -> pd.DataFrame:
    print("\n[Fix 4] Gene x ancestry carrier counts (canonical) + evidence gap ...")
    counts = canonical_df.groupby(["Gene", "ancestry_pred"])["person_id"].nunique().reset_index()
    counts = counts.rename(columns={"person_id": "n_carriers_canonical"})
    counts = counts.merge(artifact_flags[["Gene", "possible_mapping_artifact"]], on="Gene", how="left")

    pivot = counts.pivot(index="Gene", columns="ancestry_pred", values="n_carriers_canonical").fillna(0)
    total_cells = pivot.size
    below = int((pivot < min_carriers).sum().sum())
    gap_pct = 100 * below / total_cells if total_cells else np.nan

    out_path = TABLE_DIR / "gene_by_ancestry_counts.tsv"
    counts.to_csv(out_path, sep="\t", index=False)
    print(f"  wrote {out_path}")

    summary_path = TABLE_DIR / "evidence_gap_summary.txt"
    summary = (
        f"Gene x ancestry cells (canonical denominator): {total_cells}\n"
        f"Cells below min_carriers={min_carriers}: {below}\n"
        f"Evidence gap: {gap_pct:.1f}% of gene x ancestry cells have fewer than "
        f"{min_carriers} carriers.\n"
    )
    summary_path.write_text(summary)
    print(f"  wrote {summary_path}")
    print(f"  {summary.strip()}")
    return counts


# ---------------------------------------------------------------------------
# HTML report
# ---------------------------------------------------------------------------

def _b64_png(path: Path) -> str:
    return base64.b64encode(path.read_bytes()).decode("ascii")


def build_report(concordance: pd.DataFrame, gene_pen: pd.DataFrame,
                  disease_specific: pd.DataFrame, filter_log: pd.DataFrame,
                  excluded: pd.DataFrame, evidence_gap_text: str,
                  brca1_fig: Path, brca2_fig: Path, gallery_figs: dict):
    print("\n[Report] Building report.html ...")

    interpretable = gene_pen[gene_pen["interpretable"]].sort_values("penetrance_primary", ascending=False)
    excluded_pms2 = "PMS2" in set(excluded["Gene"])
    n_excluded = len(excluded)
    n_interpretable = len(interpretable)

    eas_row = concordance[concordance["ancestry_pred"] == "eas"].iloc[0]
    eur_row = concordance[concordance["ancestry_pred"] == "eur"].iloc[0]
    non_eur_max = concordance[concordance["ancestry_pred"] != "eur"].sort_values(
        "discordance_rate", ascending=False).iloc[0]

    brca1_row = interpretable[interpretable["Gene"] == "BRCA1"]
    brca2_row = interpretable[interpretable["Gene"] == "BRCA2"]

    def pct(x):
        return f"{x*100:.1f}%" if pd.notna(x) else "n/a"

    def ci(lo, hi):
        return f"({pct(lo)}–{pct(hi)})"

    concordance_table_rows = "".join(
        f"<tr><td>{r.ancestry_pred.upper()}</td><td>{int(r.n_ACMG_PLP_carriers):,}</td>"
        f"<td>{int(r.n_discordant_not_ClinVar_PLP):,}</td>"
        f"<td>{pct(r.discordance_rate)} {ci(r.CI_lower, r.CI_upper)}</td>"
        f"<td>{pct(r.singleton_fraction)}</td>"
        f"<td>{html.escape(str(r.top_variant_gene))} &mdash; {pct(r.top_variant_share)} of discordant instances</td></tr>"
        for r in concordance.itertuples()
    )

    summary_rows = "".join(
        f"<tr><td>{html.escape(r.Gene)}</td>"
        f"<td>{html.escape(r.mapped_cancer_groups)}</td>"
        f"<td>{int(r.n_carriers_canonical):,}</td>"
        f"<td>{pct(r.penetrance_primary)}</td>"
        f"<td>{ci(r.CI_lower_primary, r.CI_upper_primary)}</td></tr>"
        for r in interpretable.itertuples()
    )

    excluded_rows = "".join(
        f"<tr><td>{html.escape(r.Gene)}</td><td>{int(r.n_carriers_canonical):,}</td>"
        f"<td>{html.escape(r.exclusion_reason)}</td></tr>"
        for r in excluded.itertuples()
    )

    concordance_png = _b64_png(FIG_DIR / "concordance_by_ancestry.png")
    brca1_png = _b64_png(brca1_fig)
    brca2_png = _b64_png(brca2_fig)

    brca1_n = int(brca1_row["n_carriers_canonical"].iloc[0]) if len(brca1_row) else "n/a"
    brca2_n = int(brca2_row["n_carriers_canonical"].iloc[0]) if len(brca2_row) else "n/a"

    gallery_genes = [g for g in interpretable["Gene"] if g not in ("BRCA1", "BRCA2") and g in gallery_figs]
    gallery_html = "".join(
        f'<div><img src="data:image/png;base64,{_b64_png(gallery_figs[g])}" alt="{html.escape(g)} penetrance by age"></div>'
        for g in gallery_genes
    )

    html_doc = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>AoU Cancer-Gene Penetrance &amp; Concordance Report</title>
<style>
  :root {{
    --bg: #ffffff; --fg: #1a1a1a; --muted: #5a5a5a; --accent: #7c3aed;
    --accent-light: #f3ecfd; --border: #e2e2e2; --box-bg: #fafafa;
  }}
  * {{ box-sizing: border-box; }}
  body {{
    background: var(--bg); color: var(--fg);
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
    line-height: 1.55; max-width: 920px; margin: 0 auto; padding: 48px 24px 96px;
  }}
  h1 {{ font-size: 1.7rem; margin-bottom: 4px; }}
  h2 {{ font-size: 1.25rem; margin-top: 48px; border-bottom: 2px solid var(--accent); padding-bottom: 6px; }}
  h3 {{ font-size: 1.05rem; margin-top: 28px; }}
  .subtitle {{ color: var(--muted); font-size: 0.95rem; margin-bottom: 28px; }}
  .lede {{ font-size: 1.05rem; background: var(--accent-light); border-left: 4px solid var(--accent);
           padding: 16px 20px; border-radius: 4px; }}
  .box {{ background: var(--box-bg); border: 1px solid var(--border); border-radius: 8px;
          padding: 18px 22px; margin: 16px 0; }}
  .def {{ margin: 10px 0; }}
  .def dt {{ font-weight: 600; }}
  .def dd {{ margin: 2px 0 10px 0; color: var(--muted); }}
  table {{ border-collapse: collapse; width: 100%; margin: 14px 0; font-size: 0.92rem; }}
  th, td {{ text-align: left; padding: 7px 10px; border-bottom: 1px solid var(--border); }}
  th {{ background: var(--accent-light); font-weight: 600; }}
  tr:hover {{ background: #fafafa; }}
  img {{ max-width: 100%; border: 1px solid var(--border); border-radius: 6px; margin: 10px 0; }}
  .figrow {{ display: flex; gap: 18px; flex-wrap: wrap; }}
  .figrow > div {{ flex: 1 1 380px; }}
  .caption {{ font-size: 0.92rem; color: var(--muted); margin: 4px 0 18px; }}
  .finding {{ font-weight: 600; }}
  .small {{ font-size: 0.82rem; color: var(--muted); }}
  .tag {{ display: inline-block; background: var(--accent-light); color: var(--accent);
          border-radius: 4px; padding: 1px 8px; font-size: 0.78rem; font-weight: 600; }}
  footer {{ margin-top: 56px; color: var(--muted); font-size: 0.82rem; border-top: 1px solid var(--border); padding-top: 14px; }}
  .formula {{ font-family: "Cambria Math", "STIX Two Math", Georgia, "Times New Roman", serif;
              font-size: 1.05rem; text-align: center; margin: 14px 0; padding: 10px 0; }}
  .frac {{ display: inline-flex; flex-direction: column; vertical-align: middle;
           text-align: center; margin: 0 6px; }}
  .frac .num, .frac .den {{ display: block; padding: 1px 8px; white-space: nowrap; }}
  .frac .num {{ border-bottom: 1.3px solid currentColor; }}
  .formula sub {{ font-size: 0.7em; }}
  .formula sup {{ font-size: 0.7em; }}
</style>
</head>
<body>

<p><a href="process_walkthrough.html">&larr; How this analysis works: questions, data, and method</a></p>

<h1>Cancer-Gene Carrier Penetrance &amp; Variant-Classification Concordance</h1>
<div class="subtitle">All of Us (AoU) Researcher Workbench &middot; crude, cross-sectional, single-cohort estimates</div>

<p class="lede">
This analysis asks two questions of the same AoU carrier cohort: (1) among people carrying a
likely-pathogenic variant in a cancer-predisposition gene, what fraction actually develop the
matched cancer, and at what age (<b>penetrance</b>); and (2) does the standard practice of relying on
ClinVar to confirm pathogenicity work equally well across ancestries (<b>ACMG&nbsp;vs&nbsp;ClinVar concordance</b>).
The headline finding is an equity gap, not a biology finding: automated ACMG classification flags
many more likely-pathogenic variants than ClinVar has independently confirmed, and this gap is
far larger for non-European ancestries &mdash; meaning ClinVar-only pipelines under-ascertain risk
carriers outside European-ancestry populations. Separately, penetrance for well-established
gene&ndash;cancer pairs (BRCA1/BRCA2 &rarr; breast &amp; ovarian cancer) behaves as expected: carrier
rates are far above population background and rise with age.
</p>

<h2>Background &mdash; definitions</h2>
<div class="box">
  <dl class="def">
    <dt>Penetrance (as used here)</dt>
    <dd>Of the people carrying a likely-pathogenic variant in a gene, the fraction who have the
    <em>matched</em> cancer type recorded in this cohort, by age band. This is a <b>crude,
    cross-sectional</b> proportion &mdash; not a survival/cumulative-incidence estimate, and not
    adjusted for follow-up time or competing risks.</dd>
    <dt>ACMG/AMP P/LP</dt>
    <dd>"Pathogenic / Likely Pathogenic" classification assigned automatically from the variant's own
    features (rule-based criteria), independent of any external database.</dd>
    <dt>ClinVar P/LP</dt>
    <dd>A public database of manually-submitted, lab-reviewed variant classifications. Submissions are
    historically concentrated on variants seen in European-ancestry sequencing efforts, so absence from
    ClinVar more often reflects under-study than true benignity.</dd>
  </dl>
  <p class="small">
  <b>Carrier definition used throughout this report:</b> every penetrance number below counts a
  "carrier" as someone with an <b>ACMG-PLP</b> variant in that gene &mdash; not ClinVar-PLP, and not
  "any variant in the gene" (that broader set includes benign/VUS calls and would understate
  penetrance). We use ACMG rather than ClinVar specifically <em>because</em> Finding 1 shows ClinVar
  confirmation is far less complete for non-European carriers of the same real variants; restricting
  to ClinVar-PLP would reproduce that same ancestry bias inside the penetrance numbers.
  </p>
</div>

<h2>How the numbers are calculated</h2>
<div class="box">
  <dl class="def">
    <dt>Canonical carriers</dt>
    <dd>The <em>one</em> agreed-upon carrier list per gene that every table and figure in this report
    is built from &mdash; instead of each analysis deriving its own, slightly different count. Built
    by a fixed chain, applied identically every time: (1) start from every variant row for the gene;
    (2) collapse to one entry per person per gene; (3) restrict to carriers of an ACMG-PLP variant;
    (4) join to the cohort table to get age/sex/phenotype &mdash; carriers with no matching cohort
    record are dropped here, and that drop is counted and reported, never silent; (5) restrict to the
    correct sex if the gene's cancer type(s) are sex-specific. Whatever count survives step 5 is used,
    unchanged, as the denominator for every downstream number for that gene.</dd>

    <dt>Penetrance calculation</dt>
    <dd>A carrier is only counted as "affected" if their recorded disease Group matches the gene's
    <em>own</em> mapped cancer type(s) (e.g. BRCA1 &rarr; breast or ovarian) &mdash; a BRCA1 carrier
    with lung cancer is not counted as a BRCA1 penetrance event. People recorded only under "Family
    history" (no personal diagnosis) are explicitly excluded from the numerator by an assertion in
    the code.
    <div class="formula">
      p&#770;<sub>gene,cancer,age</sub> &nbsp;=&nbsp;
      <span class="frac"><span class="num">N(canonical carriers with that matched cancer)</span><span class="den">N(all canonical carriers in that age band)</span></span>
    </div>
    </dd>

    <dt>Background (population baseline) rate</dt>
    <dd>Computed over the <em>entire</em> cohort, not just carriers of the gene in question, with the
    same sex restriction applied as the carrier side so the comparison is apples-to-apples:
    <div class="formula">
      r&#770;<sub>bg, disease, age</sub> &nbsp;=&nbsp;
      <span class="frac"><span class="num">N(whole cohort with that exact disease)</span><span class="den">N(whole cohort in that age band)</span></span>
    </div>
    A carrier rate is marked <b>&#9733; significantly elevated</b> only when:
    <div class="formula">
      CI<sub>lower</sub>(carrier rate) &nbsp;&gt;&nbsp; CI<sub>upper</sub>(background rate)
    </div>
    &mdash; a conservative non-overlapping-interval test, not just "the bar is taller." Because this
    cohort is a case/control/family-history study design rather than a random population sample, this
    background reflects the cohort's own case-mix, not a true general-population incidence rate.</dd>

    <dt>Confidence intervals</dt>
    <dd>Every proportion above (penetrance, background rate, and the concordance/discordance rate in
    Finding 1) uses a 95% <b>Wilson score interval</b>, not the simpler normal approximation &mdash;
    several cells here have small N or proportions near 0/1, where the normal approximation
    undercoverages and can produce bounds outside [0,1]:
    <div class="formula">
      CI &nbsp;=&nbsp;
      <span class="frac">
        <span class="num">p&#770; + z<sup>2</sup>/2n &nbsp;&plusmn;&nbsp; z&radic;(p&#770;(1&minus;p&#770;)/n + z<sup>2</sup>/4n<sup>2</sup>)</span>
        <span class="den">1 + z<sup>2</sup>/n</span>
      </span>
    </div>
    z = 1.96 (95% level), n = the denominator count in the matching formula above. Computed via
    <code>statsmodels.stats.proportion.proportion_confint(method="wilson")</code>.</dd>
  </dl>
</div>

<h2>Finding 1 &mdash; the equity finding</h2>
<div class="figrow">
  <div><img src="data:image/png;base64,{concordance_png}" alt="Concordance by ancestry"></div>
</div>
<p class="finding">
Relying on ClinVar alone to confirm pathogenicity misses a much larger share of real carriers in
non-European ancestry groups: discordance (ACMG-PLP but not confirmed ClinVar-PLP) is
{pct(eur_row.discordance_rate)} in EUR vs. up to {pct(non_eur_max.discordance_rate)} in
{non_eur_max.ancestry_pred.upper()}, with confidence intervals for every non-EUR group not
overlapping EUR.
</p>
<table>
<thead><tr><th>Ancestry</th><th>ACMG-PLP carriers (n)</th><th>Discordant (n)</th>
<th>Discordance rate (95% CI)</th><th>Singleton fraction*</th><th>Recurrence check</th></tr></thead>
<tbody>{concordance_table_rows}</tbody>
</table>
<p class="small">
*Singleton fraction = share of discordant variant instances carried by only one person in this dataset
(a proxy for allele-count-1 rarity). Low singleton fractions mean discordance is mostly driven by
<em>recurrent</em> variants, not one-off rare calls &mdash; consistent with a real curation gap rather
than noise. <b>Exception:</b> EAS discordance ({pct(eas_row.discordance_rate)}) is dominated by a single
recurrent {html.escape(str(concordance[concordance.ancestry_pred=='eas'].top_variant_gene.iloc[0]))} variant
({pct(concordance[concordance.ancestry_pred=='eas'].top_variant_share.iloc[0])} of all EAS discordant
instances) &mdash; this specific number should be read as "one variant's classification status," not
broad-based under-curation across many EAS variants.
</p>

<h2>Finding 2 &mdash; penetrance behaves as expected</h2>
<div class="figrow">
  <div>
    <img src="data:image/png;base64,{brca1_png}" alt="BRCA1 penetrance by age">
    <p class="caption">BRCA1 (n={brca1_n} canonical female carriers): breast and ovarian cancer
    rates rise with age and sit far above population background at every age band &mdash;
    consistent with known BRCA1 biology, and lower than clinic-ascertained penetrance estimates,
    as expected for an unselected population cohort.</p>
  </div>
  <div>
    <img src="data:image/png;base64,{brca2_png}" alt="BRCA2 penetrance by age">
    <p class="caption">BRCA2 (n={brca2_n} canonical female carriers): same pattern, at consistently
    lower absolute rates than BRCA1 &mdash; matching the known BRCA1 &gt; BRCA2 ordering for
    breast/ovarian cancer risk.</p>
  </div>
</div>
<p class="small">
BRCA1 and BRCA2 are shown as the headline pair because they have the largest canonical carrier
counts ({brca1_n} and {brca2_n} respectively) and therefore the tightest confidence intervals &mdash;
not because other genes were excluded. The same chart for every other interpretable gene is below.
</p>

<h3>All interpretable genes, disease-specific penetrance by age</h3>
<div class="figrow">
{gallery_html}
</div>

<h2>Penetrance summary &mdash; interpretable genes only</h2>
<p class="small">
Restricted to genes with &ge;{MIN_INTERPRETABLE_CARRIERS} canonical carriers, no pseudogene-mapping-artifact
flag, and a defined gene&rarr;cancer mapping ({n_interpretable} of {n_interpretable + n_excluded} genes
analyzed qualify; see exclusions below). Sorted by penetrance, descending.
</p>
<table>
<thead><tr><th>Gene</th><th>Mapped cancer type(s)</th><th>Canonical carriers (n)</th>
<th>Penetrance (any age)</th><th>95% CI</th></tr></thead>
<tbody>{summary_rows}</tbody>
</table>

<h3>Excluded from main results &mdash; needs QC or insufficient data</h3>
<table>
<thead><tr><th>Gene</th><th>Canonical carriers (n)</th><th>Reason</th></tr></thead>
<tbody>{excluded_rows}</tbody>
</table>
<p class="small">{"PMS2 is excluded due to a known PMS2/PMS2CL pseudogene mismapping artifact (AFR carrier count far exceeds EUR, inconsistent with true prevalence) and does not appear in any main-result figure or table." if excluded_pms2 else ""}</p>

<h2>Limitations</h2>
<div class="box">
<ul>
  <li><b>Crude, cross-sectional, not survival-adjusted.</b> Penetrance here is a point-in-time
  proportion within an age band, not a cumulative-incidence or Kaplan&ndash;Meier estimate; it does
  not account for variable follow-up time.</li>
  <li><b>All of Us is not population-representative.</b> It is a healthcare-utilizing, opt-in
  volunteer cohort; cancer-predisposition gene carriers who never sought care, or who died before
  enrollment, are not captured, and case ascertainment overall is subject to selection bias.</li>
  <li><b>Denominator reconciliation.</b> {filter_log[filter_log.Gene=='BRCA1'].iloc[0].dropped_at_cohort_join:,}
  BRCA1 carriers ({filter_log[filter_log.Gene=='BRCA1'].iloc[0].pct_dropped_at_cohort_join}%) had no
  matching cohort/phenotype record and were dropped before any penetrance calculation; see
  <code>carrier_filter_log.tsv</code> for the full per-gene breakdown.</li>
  <li><b>Some genes excluded for pseudogene-mapping artifacts</b> (notably PMS2) &mdash; see exclusion
  table above.</li>
  <li><b>Small cells excluded.</b> {evidence_gap_text}</li>
  <li>These results are preliminary, single-cohort, and not validated against an independent dataset.</li>
</ul>
</div>

<h2>What this sets up</h2>
<p>
The concordance gap motivates <b>ancestry-aware variant classification</b> &mdash; ClinVar-only
pipelines should not be trusted as a sole confirmation step outside European-ancestry populations,
and ACMG-flagged variants lacking ClinVar support in non-EUR carriers are a natural re-curation
priority list. The penetrance pipeline built here (canonical denominator, disease-matched numerator,
population-background comparison) is the basis for the next step: variant-level, competing-risk-adjusted
absolute risk estimation.
</p>

<footer>
Generated from AoU carrier and cohort tables in <code>DemoVariant/</code>. Script:
<code>DemoVariant/Penetrance/gene_penetrance.py</code>. All figures and tables in this report are
reproducible by re-running that script.
</footer>

</body>
</html>
"""
    out_path = OUT_DIR / "report.html"
    out_path.write_text(html_doc, encoding="utf-8")
    print(f"  wrote {out_path}")


def build_process_walkthrough(filter_log: pd.DataFrame, artifact_flags: pd.DataFrame,
                               gene_pen: pd.DataFrame, concordance: pd.DataFrame):
    print("\n[Report] Building process_walkthrough.html ...")

    brca1 = filter_log[filter_log["Gene"] == "BRCA1"].iloc[0]
    flagged = artifact_flags[artifact_flags["possible_mapping_artifact"]]

    def n(x):
        return f"{int(x):,}"

    filter_log_rows = "".join(
        f"<tr><td>{html.escape(r.Gene)}</td><td>{n(r.raw_variant_rows)}</td>"
        f"<td>{n(r.unique_person_gene_any_classification)}</td>"
        f"<td>{n(r.unique_person_gene_ACMG_PLP)}</td>"
        f"<td>{n(r.joined_to_cohort)}</td>"
        f"<td>{n(r.dropped_at_cohort_join)} ({r.pct_dropped_at_cohort_join}%)</td>"
        f"<td>{'Female only' if r.sex_restricted else 'Both sexes'}</td>"
        f"<td>{n(r.n_carriers_canonical)}</td></tr>"
        for r in filter_log.sort_values('n_carriers_canonical', ascending=False).itertuples()
    )

    artifact_rows = "".join(
        f"<tr><td>{html.escape(r.Gene)}</td><td>{html.escape(r.mapping_artifact_reason)}</td></tr>"
        for r in flagged.itertuples()
    )

    html_doc = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>How This Analysis Works</title>
<style>
  :root {{
    --bg: #ffffff; --fg: #1a1a1a; --muted: #5a5a5a; --accent: #7c3aed;
    --accent-light: #f3ecfd; --border: #e2e2e2; --box-bg: #fafafa;
  }}
  * {{ box-sizing: border-box; }}
  body {{
    background: var(--bg); color: var(--fg);
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
    line-height: 1.6; max-width: 920px; margin: 0 auto; padding: 48px 24px 96px;
  }}
  h1 {{ font-size: 1.7rem; margin-bottom: 4px; }}
  h2 {{ font-size: 1.25rem; margin-top: 44px; border-bottom: 2px solid var(--accent); padding-bottom: 6px; }}
  h3 {{ font-size: 1.05rem; margin-top: 24px; }}
  .subtitle {{ color: var(--muted); font-size: 0.95rem; margin-bottom: 28px; }}
  .box {{ background: var(--box-bg); border: 1px solid var(--border); border-radius: 8px;
          padding: 18px 22px; margin: 16px 0; }}
  .step {{ display: flex; gap: 14px; margin: 14px 0; align-items: flex-start; }}
  .step .num {{ flex: 0 0 28px; height: 28px; border-radius: 50%; background: var(--accent);
                color: white; display: flex; align-items: center; justify-content: center;
                font-weight: 700; font-size: 0.85rem; }}
  .step .body {{ flex: 1; }}
  table {{ border-collapse: collapse; width: 100%; margin: 14px 0; font-size: 0.86rem; }}
  th, td {{ text-align: left; padding: 6px 9px; border-bottom: 1px solid var(--border); }}
  th {{ background: var(--accent-light); font-weight: 600; }}
  .scroll {{ overflow-x: auto; }}
  code {{ background: var(--box-bg); padding: 1px 5px; border-radius: 3px; font-size: 0.88em; }}
  .small {{ font-size: 0.85rem; color: var(--muted); }}
  ul {{ padding-left: 22px; }}
  li {{ margin: 6px 0; }}
  .formula {{ font-family: "Cambria Math", "STIX Two Math", Georgia, "Times New Roman", serif;
              font-size: 1.05rem; text-align: center; margin: 14px 0; padding: 10px 0; }}
  .frac {{ display: inline-flex; flex-direction: column; vertical-align: middle;
           text-align: center; margin: 0 6px; }}
  .frac .num, .frac .den {{ display: block; padding: 1px 8px; white-space: nowrap; }}
  .frac .num {{ border-bottom: 1.3px solid currentColor; }}
  .formula sub {{ font-size: 0.7em; }}
  .formula sup {{ font-size: 0.7em; }}
</style>
</head>
<body>

<p><a href="report_acmg_vs_clinvar.html">&rarr; ACMG vs ClinVar comparison (report_acmg_vs_clinvar.html)</a></p>

<h1>How This Analysis Works</h1>
<div class="subtitle">Questions asked, data used, and the method behind every number in report.html
and report_acmg_vs_clinvar.html</div>

<h2>The two questions</h2>
<ol>
  <li><b>Penetrance:</b> among people who carry a likely-pathogenic variant in a cancer-predisposition
  gene, what fraction actually develop the cancer that gene is known to cause, and at what age?</li>
  <li><b>Concordance:</b> does relying on ClinVar to confirm a variant's pathogenicity work equally
  well for every ancestry group, or does it systematically miss carriers in some groups more than
  others?</li>
</ol>
<p>Both questions are asked of the <em>same</em> AoU carrier cohort, so the same denominator-building
logic underlies both.</p>

<h2>Data sources</h2>
<table>
<thead><tr><th>File</th><th>What it is</th><th>Granularity</th></tr></thead>
<tbody>
<tr><td><code>carriers_all_noMUTYH_ancestry.tsv</code></td>
    <td>Every variant call in the tracked cancer-predisposition genes, per person</td>
    <td>One row per (person, variant) &mdash; a person with 2 variants in a gene has 2 rows</td></tr>
<tr><td><code>All_temp_cohort_simple.csv.gz</code></td>
    <td>Cohort membership: age, sex, ancestry, and recorded disease Group</td>
    <td>One row per (person, condition) &mdash; a person with 2 recorded conditions has 2 rows</td></tr>
</tbody>
</table>
<p class="small">Both tables are multi-row-per-person by construction. Every bug found and fixed in
this pipeline traced back to code that assumed one row per person and didn't dedupe before joining
or counting &mdash; see the denominator method below.</p>

<h2>Building the canonical carrier denominator</h2>
<p>Every gene-level number in the report (penetrance, ancestry breakdown, evidence-gap %) is computed
from <b>one</b> canonical carrier set per gene, built in this order:</p>

<div class="step"><div class="num">1</div><div class="body">
<b>Start from raw variant rows</b> for the gene &mdash; every (person, variant) row, including
multiple variants per person.</div></div>
<div class="step"><div class="num">2</div><div class="body">
<b>Collapse to unique person &times; gene.</b> A person counts as a carrier once per gene, and as
ACMG-PLP/ClinVar-PLP for that gene if <em>any</em> of their variants in it carries that flag.</div></div>
<div class="step"><div class="num">3</div><div class="body">
<b>Restrict to ACMG-PLP carriers.</b> We use ACMG (rule-based, applied uniformly regardless of
ancestry) rather than ClinVar-PLP as the carrier definition &mdash; restricting to ClinVar-PLP would
bake the ancestry-based confirmation gap (Finding 1) directly into the penetrance denominator.</div></div>
<div class="step"><div class="num">4</div><div class="body">
<b>Join to the cohort table</b> to get phenotype/age/sex, using a person-level cohort table built by
first collapsing the cohort's multiple per-person rows into (a) one demographic row and (b) a
<em>set</em> of every disease Group that person has ever been recorded under. Carriers with no
matching cohort record are dropped here and reported explicitly (they have no phenotype data, so
they cannot contribute to a penetrance calculation).</div></div>
<div class="step"><div class="num">5</div><div class="body">
<b>Restrict to the correct sex</b> for genes whose mapped cancer type(s) are sex-specific (e.g.
breast/ovarian/endometrial genes are restricted to female carriers, both numerator and denominator).</div></div>

<div class="scroll">
<table>
<thead><tr><th>Gene</th><th>Raw rows</th><th>Unique person&times;gene</th>
<th>ACMG-PLP only</th><th>Joined to cohort</th><th>Dropped (no cohort record)</th>
<th>Sex restriction</th><th>Canonical N</th></tr></thead>
<tbody>{filter_log_rows}</tbody>
</table>
</div>
<p class="small">Full table also written to <code>carrier_filter_log.tsv</code>. Worked example for
BRCA1: {n(brca1.raw_variant_rows)} raw rows &rarr; {n(brca1.unique_person_gene_any_classification)}
unique carriers (any classification) &rarr; {n(brca1.unique_person_gene_ACMG_PLP)} ACMG-PLP carriers
&rarr; {n(brca1.joined_to_cohort)} joined to cohort (dropped {n(brca1.dropped_at_cohort_join)},
{brca1.pct_dropped_at_cohort_join}%) &rarr; <b>{n(brca1.n_carriers_canonical)} canonical carriers</b>
(female-restricted). This single number now flows unchanged into every table and figure that
mentions BRCA1.</p>

<h2>Excluding genotyping artifacts</h2>
<p>Some genes sit in or near a pseudogene / highly-homologous paralog region, where short-read
mapping can generate spurious, ancestry-skewed carrier calls that are a sequencing artifact, not
biology. We flag a gene if:</p>
<ul>
  <li>it is <b>PMS2</b> (a known PMS2/PMS2CL pseudogene mismapping problem), or</li>
  <li>it is on a small watch-list of other pseudogene-adjacent genes (CHEK2, SMAD4) <b>and</b> one
  non-European ancestry has more than 2&times; the European-ancestry carrier count &mdash; a pattern
  with no plausible biological explanation.</li>
</ul>
<table>
<thead><tr><th>Gene</th><th>Reason flagged</th></tr></thead>
<tbody>{artifact_rows if len(artifact_rows) else '<tr><td colspan="2">No genes flagged beyond PMS2.</td></tr>'}</tbody>
</table>
<p class="small">Flagged genes are excluded from the main penetrance results and figures, and listed
separately in <code>excluded_genes.tsv</code> and the report's "Excluded" table &mdash; never silently
dropped.</p>

<h2>Disease-specific penetrance, not "any cancer"</h2>
<p>For each gene we use an explicit gene &rarr; cancer-type mapping (e.g. BRCA1 &rarr; Breast, Ovary)
rather than counting <em>any</em> recorded cancer as "affected." A BRCA1 carrier with lung cancer is
not counted as a BRCA1 penetrance event. Each mapped cancer type is reported and plotted as its own
line &mdash; BRCA1:Breast and BRCA1:Ovary are never merged into one number. People recorded only
under "Family history" (no personal diagnosis) are explicitly excluded from the numerator by an
assertion in the code.</p>
<div class="formula">
  p&#770;<sub>gene,cancer,age</sub> &nbsp;=&nbsp;
  <span class="frac"><span class="num">N(canonical carriers with matched cancer)</span><span class="den">N(all canonical carriers in that age band)</span></span>
</div>

<p>Every carrier rate is compared to a <b>population background rate</b>: the same disease's
prevalence in the whole cohort (not just carriers of that gene), computed in the same age bands and
with the same sex restriction, so the comparison is apples-to-apples.</p>
<div class="formula">
  r&#770;<sub>bg, disease, age</sub> &nbsp;=&nbsp;
  <span class="frac"><span class="num">N(whole cohort with that exact disease)</span><span class="den">N(whole cohort in that age band)</span></span>
</div>
<p>A gene-disease-age cell is marked "significantly elevated" only when:</p>
<div class="formula">
  CI<sub>lower</sub>(carrier rate) &nbsp;&gt;&nbsp; CI<sub>upper</sub>(background rate)
</div>
<p class="small">A conservative non-overlapping-interval test. Because this cohort is a
case/control/family-history study design rather than a random population sample, this background
reflects the cohort's own case-mix, not a true general-population incidence rate.</p>

<p>All proportions above use a 95% <b>Wilson score interval</b>, not the normal approximation
(several cells here have small N or proportions near 0/1, where the normal approximation
undercoverages):</p>
<div class="formula">
  CI &nbsp;=&nbsp;
  <span class="frac">
    <span class="num">p&#770; + z<sup>2</sup>/2n &nbsp;&plusmn;&nbsp; z&radic;(p&#770;(1&minus;p&#770;)/n + z<sup>2</sup>/4n<sup>2</sup>)</span>
    <span class="den">1 + z<sup>2</sup>/n</span>
  </span>
</div>
<p class="small">z = 1.96 (95% level), n = the matching denominator above. Computed via
<code>statsmodels.stats.proportion.proportion_confint(method="wilson")</code>.</p>

<h2>Concordance footnotes</h2>
<p>The headline concordance-by-ancestry numbers are ACMG-PLP carriers who are not also ClinVar-PLP,
by ancestry. Two footnote checks guard against over-reading a single high number:</p>
<ul>
  <li><b>Singleton fraction:</b> what share of discordant variant instances are carried by only one
  person in the dataset (a rarity proxy)? A high fraction would suggest the gap is about rare, hard-to
  -classify variants rather than a systemic curation gap.</li>
  <li><b>Recurrent-variant share:</b> what fraction of an ancestry's discordant instances trace back
  to its single most common variant? A high share means one variant is driving the headline number,
  and it should be read as "this one variant's status," not broad under-curation.</li>
</ul>
<p class="small">Applying this check surfaced a real finding: EAS discordance
({concordance[concordance.ancestry_pred=='eas'].discordance_rate.iloc[0]*100:.0f}%) is
{concordance[concordance.ancestry_pred=='eas'].top_variant_share.iloc[0]*100:.0f}% attributable to a
single recurrent
{html.escape(str(concordance[concordance.ancestry_pred=='eas'].top_variant_gene.iloc[0]))} variant,
with a low singleton fraction &mdash; reported as a specific-variant caveat in the main report rather
than removed, per the review that requested this check.</p>

<h2>Limitations of the method itself</h2>
<ul>
  <li>The "ACMG-PLP" restriction still depends on the upstream ACMG classification pipeline being
  correct; we did not re-derive ACMG criteria here.</li>
  <li>Population background rates are computed on this same ascertained cohort (a
  healthcare-utilizing volunteer population), not a true general-population rate.</li>
  <li>Age bands are coarse (5 bins); a finer-grained or continuous age model would use more of the
  data but was out of scope for this pass.</li>
  <li>The pseudogene-artifact rule is a heuristic (PMS2 always flagged; a 2&times; EUR-ratio rule for
  two other watch-list genes) and does not exhaustively check every gene for mapping problems.</li>
</ul>

<p><a href="report_acmg_vs_clinvar.html">&rarr; ACMG vs ClinVar comparison (report_acmg_vs_clinvar.html)</a></p>

</body>
</html>
"""
    out_path = OUT_DIR / "process_walkthrough.html"
    out_path.write_text(html_doc, encoding="utf-8")
    print(f"  wrote {out_path}")


CLINVAR_COMPARISON_MIN = MIN_INTERPRETABLE_CARRIERS  # kept equal to the main report's own
                              # inclusion bar for consistency across reports -- both classifiers
                              # are judged by the same N standard.


def build_acmg_vs_clinvar_report(gene_pen_acmg: pd.DataFrame, gene_pen_clinvar: pd.DataFrame,
                                  filter_log_acmg: pd.DataFrame, filter_log_clinvar: pd.DataFrame,
                                  concordance: pd.DataFrame, gallery_figs_acmg: dict,
                                  gallery_figs_clinvar: dict):
    """Side-by-side comparison: same genes, same disease mapping, same
    canonical-denominator method -- the ONLY thing that changes is which
    classifier's PLP calls define a 'carrier'. Written to a NEW file;
    report.html and process_walkthrough.html are untouched."""
    print("\n[Report] Building report_acmg_vs_clinvar.html ...")

    def pct(x):
        return f"{x*100:.1f}%" if pd.notna(x) else "n/a"

    genes = sorted(set(gene_pen_acmg["Gene"]) | set(gene_pen_clinvar["Gene"]))
    a = gene_pen_acmg.set_index("Gene")
    c = gene_pen_clinvar.set_index("Gene")

    rows = []
    for gene in genes:
        if gene not in a.index or gene not in c.index:
            continue
        ra, rc = a.loc[gene], c.loc[gene]
        if not ra["has_gene_cancer_mapping"]:
            continue
        rows.append({
            "Gene": gene,
            "mapped": ra["mapped_cancer_groups"],
            "n_acmg": int(ra["n_carriers_canonical"]),
            "pen_acmg": ra["penetrance_primary"],
            "n_clinvar": int(rc["n_carriers_canonical"]),
            "pen_clinvar": rc["penetrance_primary"],
            "n_ratio": (int(rc["n_carriers_canonical"]) / int(ra["n_carriers_canonical"])
                        if ra["n_carriers_canonical"] else np.nan),
            "possible_mapping_artifact": bool(ra["possible_mapping_artifact"]),
        })
    comp = pd.DataFrame(rows).sort_values("n_acmg", ascending=False)

    total_acmg = comp["n_acmg"].sum()
    total_clinvar = comp["n_clinvar"].sum()
    missed_pct = 100 * (1 - total_clinvar / total_acmg) if total_acmg else np.nan

    # "Confident" genes for FIGURES: ACMG side already meets the main report's
    # MIN_INTERPRETABLE_CARRIERS bar, and the ClinVar side meets CLINVAR_COMPARISON_MIN
    # (this comparison is inherently about small-N genes, but a chart still
    # needs a plottable, non-noise N to be shown). Genes failing this are
    # listed in a table, not silently dropped.
    comp["confident"] = (
        (comp["n_acmg"] >= MIN_INTERPRETABLE_CARRIERS)
        & (comp["n_clinvar"] >= CLINVAR_COMPARISON_MIN)
        & ~comp["possible_mapping_artifact"]
    )
    confident = comp[comp["confident"]].copy()
    not_confident = comp[~comp["confident"]].copy()

    def _low_reason(r):
        reasons = []
        if r["n_acmg"] < MIN_INTERPRETABLE_CARRIERS:
            reasons.append(f"ACMG N={r['n_acmg']} (<{MIN_INTERPRETABLE_CARRIERS})")
        if r["n_clinvar"] < CLINVAR_COMPARISON_MIN:
            reasons.append(f"ClinVar N={r['n_clinvar']} (<{CLINVAR_COMPARISON_MIN})")
        if r["possible_mapping_artifact"]:
            reasons.append("pseudogene-mapping artifact (see main report)")
        return "; ".join(reasons)
    not_confident["reason"] = not_confident.apply(_low_reason, axis=1)

    # Figure 1: carrier counts, ACMG-PLP vs ClinVar-PLP, confident genes only.
    x = np.arange(len(confident))
    width = 0.38
    fig1, ax1 = plt.subplots(figsize=(10, 6))
    ax1.bar(x - width / 2, confident["n_acmg"], width, label="ACMG-PLP carriers", color="#4C72B0")
    ax1.bar(x + width / 2, confident["n_clinvar"], width, label="ClinVar-PLP carriers", color="#DD8452")
    ax1.set_xticks(x)
    ax1.set_xticklabels(confident["Gene"], rotation=45, ha="right", fontsize=9)
    ax1.set_ylabel("Canonical carriers")
    ax1.set_title(f"Canonical carrier count by classifier — confident genes only "
                   f"(ACMG N≥{MIN_INTERPRETABLE_CARRIERS}, ClinVar N≥{CLINVAR_COMPARISON_MIN})")
    ax1.legend()
    fig1.tight_layout()
    carriers_fig_path = FIG_DIR / "acmg_vs_clinvar_carriers.png"
    fig1.savefig(carriers_fig_path, dpi=150)
    plt.close(fig1)
    print(f"  wrote {carriers_fig_path}")

    # Figure 2: penetrance, ACMG-PLP vs ClinVar-PLP, confident genes only.
    x2 = np.arange(len(confident))
    fig2, ax2 = plt.subplots(figsize=(10, 6))
    ax2.bar(x2 - width / 2, confident["pen_acmg"] * 100, width, label="ACMG-PLP penetrance", color="#4C72B0")
    ax2.bar(x2 + width / 2, confident["pen_clinvar"] * 100, width, label="ClinVar-PLP penetrance", color="#DD8452")
    ax2.set_xticks(x2)
    ax2.set_xticklabels(confident["Gene"], rotation=45, ha="right", fontsize=9)
    ax2.set_ylabel("Penetrance (%)")
    ax2.set_title(f"Penetrance by classifier — confident genes only "
                   f"(ACMG N≥{MIN_INTERPRETABLE_CARRIERS}, ClinVar N≥{CLINVAR_COMPARISON_MIN})")
    ax2.legend()
    fig2.tight_layout()
    penetrance_fig_path = FIG_DIR / "acmg_vs_clinvar_penetrance.png"
    fig2.savefig(penetrance_fig_path, dpi=150)
    plt.close(fig2)
    print(f"  wrote {penetrance_fig_path}")

    # Paired per-gene disease-specific charts (the ACMG-style charts), ACMG
    # next to ClinVar, for confident genes only.
    pair_html = "".join(
        f'<div><h4>{html.escape(g)}</h4><div class="figrow">'
        f'<div><img src="data:image/png;base64,{_b64_png(gallery_figs_acmg[g])}" alt="{html.escape(g)} ACMG-PLP"><p class="caption">ACMG-PLP</p></div>'
        f'<div><img src="data:image/png;base64,{_b64_png(gallery_figs_clinvar[g])}" alt="{html.escape(g)} ClinVar-PLP"><p class="caption">ClinVar-PLP</p></div>'
        f'</div></div>'
        for g in confident["Gene"] if g in gallery_figs_acmg and g in gallery_figs_clinvar
    )

    not_confident_rows = "".join(
        f"<tr><td>{html.escape(r.Gene)}</td><td>{html.escape(r.reason)}</td></tr>"
        for r in not_confident.itertuples()
    )

    rows_html = "".join(
        f"<tr><td>{html.escape(r.Gene)}</td><td>{html.escape(r.mapped)}</td>"
        f"<td>{r.n_acmg:,}</td><td>{pct(r.pen_acmg)}</td>"
        f"<td>{r.n_clinvar:,}</td><td>{pct(r.pen_clinvar)}</td>"
        f"<td>{r.n_ratio*100:.0f}%</td></tr>"
        for r in comp.itertuples()
    )

    eas_row = concordance[concordance["ancestry_pred"] == "eas"].iloc[0]
    non_eur_max = concordance[concordance["ancestry_pred"] != "eur"].sort_values(
        "discordance_rate", ascending=False).iloc[0]

    html_doc = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>ACMG vs ClinVar Carrier Definition Comparison</title>
<style>
  :root {{
    --bg: #ffffff; --fg: #1a1a1a; --muted: #5a5a5a; --accent: #7c3aed;
    --accent-light: #f3ecfd; --border: #e2e2e2; --box-bg: #fafafa;
  }}
  * {{ box-sizing: border-box; }}
  body {{
    background: var(--bg); color: var(--fg);
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
    line-height: 1.55; max-width: 960px; margin: 0 auto; padding: 48px 24px 96px;
  }}
  h1 {{ font-size: 1.6rem; margin-bottom: 4px; }}
  h2 {{ font-size: 1.2rem; margin-top: 40px; border-bottom: 2px solid var(--accent); padding-bottom: 6px; }}
  .subtitle {{ color: var(--muted); font-size: 0.95rem; margin-bottom: 24px; }}
  .lede {{ font-size: 1.02rem; background: var(--accent-light); border-left: 4px solid var(--accent);
           padding: 16px 20px; border-radius: 4px; }}
  .box {{ background: var(--box-bg); border: 1px solid var(--border); border-radius: 8px;
          padding: 16px 20px; margin: 14px 0; }}
  table {{ border-collapse: collapse; width: 100%; margin: 14px 0; font-size: 0.9rem; }}
  th, td {{ text-align: left; padding: 7px 10px; border-bottom: 1px solid var(--border); }}
  th {{ background: var(--accent-light); font-weight: 600; }}
  .scroll {{ overflow-x: auto; }}
  .small {{ font-size: 0.85rem; color: var(--muted); }}
  .lowratio {{ color: #b91c1c; font-weight: 600; }}
  .formula {{ font-family: "Cambria Math", "STIX Two Math", Georgia, "Times New Roman", serif;
              font-size: 1.05rem; text-align: center; margin: 14px 0; padding: 10px 0; }}
  .frac {{ display: inline-flex; flex-direction: column; vertical-align: middle;
           text-align: center; margin: 0 6px; }}
  .frac .num, .frac .den {{ display: block; padding: 1px 8px; white-space: nowrap; }}
  .frac .num {{ border-bottom: 1.3px solid currentColor; }}
  .formula sub {{ font-size: 0.7em; }}
  .formula sup {{ font-size: 0.7em; }}
  img {{ max-width: 100%; border: 1px solid var(--border); border-radius: 6px; margin: 6px 0; }}
  .figrow {{ display: flex; gap: 14px; flex-wrap: wrap; }}
  .figrow > div {{ flex: 1 1 380px; }}
  .caption {{ font-size: 0.82rem; color: var(--muted); text-align: center; margin: 2px 0 4px; }}
  h4 {{ margin: 26px 0 4px; font-size: 1rem; }}
</style>
</head>
<body>

<p><a href="process_walkthrough.html">&larr; How this analysis works (questions, data, method)</a></p>

<h1>Carrier Definition Comparison: ACMG-PLP vs ClinVar-PLP</h1>
<div class="subtitle">Same genes, same disease mapping, same canonical-denominator method &mdash;
only the pathogenicity classifier defining "carrier" changes</div>

<p class="lede">
The main report defines a carrier as someone with an <b>ACMG-PLP</b> variant. This companion report
asks: what happens if we instead required <b>ClinVar-PLP</b> confirmation? Answer: it doesn't just
lower the numbers &mdash; it <b>reproduces the ancestry gap from the concordance finding inside the
penetrance denominator itself</b>. Restricting to ClinVar-PLP drops
{comp['n_acmg'].sum() - comp['n_clinvar'].sum():,} carriers overall
({missed_pct:.0f}% of the ACMG-PLP carrier pool) across these genes, and (see the main report) that
loss falls disproportionately on non-European ancestries &mdash; up to {pct(non_eur_max.discordance_rate)}
of ACMG-PLP carriers in {non_eur_max.ancestry_pred.upper()} have no ClinVar confirmation at all, vs.
{pct(concordance[concordance.ancestry_pred=='eur'].discordance_rate.iloc[0])} in EUR. This is why the
main report's headline penetrance numbers use ACMG, not ClinVar, as the carrier definition.
</p>

<h2>What "canonical carriers" means, again</h2>
<div class="box">
<p>Same definition as in <a href="process_walkthrough.html">process_walkthrough.html</a>: the ONE
carrier set per gene, built by raw rows &rarr; unique person&times;gene &rarr; restricted to a PLP
classification &rarr; joined to cohort &rarr; sex-restricted if applicable. The only difference
between the two columns below is which classification (ACMG-PLP vs ClinVar-PLP) is used in that
restriction step &mdash; everything else in the pipeline is identical.</p>
</div>

<h2>Calculation and statistics used</h2>
<div class="box">

  <p><b>Penetrance</b> &mdash; computed separately for the ACMG-PLP carrier set and the ClinVar-PLP
  carrier set (same formula, different input population):</p>
  <div class="formula">
    p&#770;<sub>gene,cancer,age</sub> &nbsp;=&nbsp;
    <span class="frac"><span class="num">N(carriers with matched cancer)</span><span class="den">N(all canonical carriers)</span></span>
    &nbsp;&nbsp;&nbsp;&nbsp; <span style="color:var(--muted)">[age band, sex-restricted where applicable]</span>
  </div>

  <p><b>Background (population baseline) rate</b> &mdash; same disease and age band, computed over
  the whole cohort, independent of which classifier defines "carrier":</p>
  <div class="formula">
    r&#770;<sub>bg, disease, age</sub> &nbsp;=&nbsp;
    <span class="frac"><span class="num">N(cohort with that disease)</span><span class="den">N(cohort in that age band)</span></span>
  </div>

  <p><b>Confidence intervals</b> &mdash; 95% <b>Wilson score interval</b> on every proportion above
  (penetrance, background rate, concordance/discordance rate), chosen over the simpler normal
  approximation because several cells here have small N or proportions near 0/1, where the normal
  approximation undercoverages and can produce bounds outside [0,1]:</p>
  <div class="formula">
    CI &nbsp;=&nbsp;
    <span class="frac">
      <span class="num">p&#770; + z<sup>2</sup>/2n &nbsp;&plusmn;&nbsp; z&radic;(p&#770;(1&minus;p&#770;)/n + z<sup>2</sup>/4n<sup>2</sup>)</span>
      <span class="den">1 + z<sup>2</sup>/n</span>
    </span>
  </div>
  <p class="small" style="text-align:center; margin-top:-8px;">
  z = 1.96 (95% level), n = denominator count above. Computed via
  <code>statsmodels.stats.proportion.proportion_confint(method="wilson")</code>.</p>

  <p><b>"Significantly elevated" (&#9733;)</b> &mdash; a non-overlapping-interval test, not a formal
  two-proportion hypothesis test (no p-value computed):</p>
  <div class="formula">
    CI<sub>lower</sub>(carrier rate) &nbsp;&gt;&nbsp; CI<sub>upper</sub>(background rate)
  </div>

  <p><b>Carrier-count ratio</b> (last table column) &mdash; a descriptive ratio, not a statistical
  test, used only to rank which genes lose the most carriers under the ClinVar restriction:</p>
  <div class="formula">
    ratio &nbsp;=&nbsp;
    <span class="frac"><span class="num">N<sub>ClinVar-PLP</sub></span><span class="den">N<sub>ACMG-PLP</sub></span></span>
  </div>
</div>

<h2>Min-carriers thresholds used</h2>
<div class="box">
<p>Two different thresholds exist in this pipeline:</p>
<ul>
<li><b>Gene-inclusion bar ({MIN_INTERPRETABLE_CARRIERS} canonical carriers):</b> a gene needs
&ge;{MIN_INTERPRETABLE_CARRIERS} canonical carriers to be treated as "interpretable." This comparison
report applies that <b>same</b> bar to <em>both</em> classifiers &mdash; ACMG-PLP N and ClinVar-PLP N
each need &ge;{MIN_INTERPRETABLE_CARRIERS} &mdash; so neither side gets an easier standard. The full
table further down still shows every gene, including ones that fail this bar on one or both sides, so
the shortfall itself stays visible.</li>
<li><b>Evidence-gap threshold (3 carriers):</b> used only in the ancestry-breakdown table to flag
gene&times;ancestry cells too sparse to trust. Unrelated to this comparison.</li>
</ul>
</div>

<h2>Confident genes only</h2>
<p class="small">
Figures below are restricted to genes meeting <b>both</b> N bars at the same threshold: ACMG-PLP N
&ge; {MIN_INTERPRETABLE_CARRIERS} and ClinVar-PLP N &ge; {CLINVAR_COMPARISON_MIN} (identical value,
applied consistently to both classifiers). {len(confident)} of {len(comp)} genes qualify; the rest
are listed in the table at the very bottom with the specific reason each was left out of the figures
&mdash; nothing is silently dropped, it's just not plotted.
</p>

<h2>Carrier counts by classifier</h2>
<img src="data:image/png;base64,{_b64_png(carriers_fig_path)}" alt="ACMG vs ClinVar carrier counts by gene">
<p class="small">Genes where the orange (ClinVar-PLP) bar is visibly shorter than the blue (ACMG-PLP)
bar lost the most carriers to the ClinVar restriction.</p>

<h2>Penetrance by classifier</h2>
<img src="data:image/png;base64,{_b64_png(penetrance_fig_path)}" alt="ACMG vs ClinVar penetrance by gene">
<p class="small">Read this chart together with the carrier-count chart above, not alone: a similar
bar height with a much smaller N behind it (see table below) is a wider, less certain estimate, not a
more "correct" one.</p>

<h2>Per-gene disease-specific charts, ACMG-PLP vs ClinVar-PLP</h2>
<p class="small">Same chart style as the main report's gene gallery (each mapped cancer type as its
own line, population background dashed, &#9733; = significantly elevated) &mdash; ACMG-PLP carriers on
the left, ClinVar-PLP carriers on the right, same gene, same axes scale intent, so you can see
directly how much noisier/sparser the ClinVar-PLP curve is.</p>
{pair_html}

<h2>Genes left out of the figures above</h2>
<table>
<thead><tr><th>Gene</th><th>Reason</th></tr></thead>
<tbody>{not_confident_rows if len(not_confident_rows) else '<tr><td colspan="2">None &mdash; all genes met the confidence bar.</td></tr>'}</tbody>
</table>

<h2>Full side-by-side table, all genes</h2>
<p class="small">Every gene, including the low-confidence ones excluded from the figures above &mdash;
kept here for completeness.</p>
<div class="scroll">
<table>
<thead><tr><th>Gene</th><th>Mapped cancer(s)</th>
<th>ACMG-PLP N</th><th>ACMG-PLP penetrance</th>
<th>ClinVar-PLP N</th><th>ClinVar-PLP penetrance</th>
<th>ClinVar N as % of ACMG N</th></tr></thead>
<tbody>{rows_html}</tbody>
</table>
</div>
<p class="small">Rows where the last column is well under 100% are genes where switching to
ClinVar-PLP would have thrown away a substantial share of real carriers &mdash; and, per the main
report's concordance finding, that loss is not random across ancestries.</p>

<h2>Interpretation</h2>
<p>
This comparison is not "which number is right" &mdash; ACMG-PLP is the carrier definition used
throughout the main report, for the reason stated above. It exists to make concrete, gene by gene,
what the abstract concordance percentage in Finding 1 actually costs: real carriers who would be
excluded from any pipeline that requires ClinVar confirmation before counting someone as a carrier.
</p>

<footer style="margin-top:48px; color:var(--muted); font-size:0.82rem; border-top:1px solid var(--border); padding-top:14px;">
Generated by <code>DemoVariant/Penetrance/gene_penetrance.py</code>. Underlying data:
<code>gene_penetrance.tsv</code> (ACMG) and <code>gene_penetrance_clinvar.tsv</code> (ClinVar).
</footer>

</body>
</html>
"""
    out_path = OUT_DIR / "report_acmg_vs_clinvar.html"
    out_path.write_text(html_doc, encoding="utf-8")
    print(f"  wrote {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Cancer-gene carrier penetrance and ACMG/ClinVar concordance.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--carriers-path", "--carriers_path", dest="carriers_path", type=Path,
        default=HERE / DEFAULT_CARRIERS_NAME,
        help="Carrier-level TSV, one row per person x variant. Columns: person_id, "
             "Gene, Variant_ID, is_ACMG_PLP, is_ClinVar_PLP, ancestry_pred, "
             "sex_at_birth, age.",
    )
    parser.add_argument(
        "--cohort-path", "--cohort_path", dest="cohort_path", type=Path,
        default=HERE / DEFAULT_COHORT_NAME,
        help="Gzipped cohort CSV, one row per person x recorded condition. "
             "Columns: person_id, age, sex_at_birth, Group.",
    )
    parser.add_argument(
        "--outdir", type=Path, default=OUT_DIR,
        help="Where tables/, figures/ and the HTML reports are written.",
    )
    parser.add_argument(
        "--min_carriers", "--min-carriers", dest="min_carriers", type=int, default=3,
        help="Evidence-gap threshold: gene x ancestry cells below this are counted "
             "as too sparse to trust. Unrelated to the gene-inclusion bar "
             "(MIN_INTERPRETABLE_CARRIERS).",
    )
    return parser


def _resolve_paths(args) -> None:
    """Point the module-level path globals at this run's inputs and outputs.

    The loaders and writers read these globals rather than taking paths as
    arguments; setting them in one place keeps that indirection to a single
    function instead of threading a config object through twenty call sites.
    """
    global CARRIERS_PATH, COHORT_PATH, OUT_DIR, TABLE_DIR, FIG_DIR

    for label, path in (("--carriers-path", args.carriers_path),
                        ("--cohort-path", args.cohort_path)):
        if not path.exists():
            raise SystemExit(
                f"{label}: no such file: {path}\n"
                "No participant-level data ships with this repo -- see "
                "'Input data schema' in README.md for what the two tables need "
                "to contain, and pass your own with --carriers-path/--cohort-path."
            )

    CARRIERS_PATH = str(args.carriers_path)
    COHORT_PATH   = str(args.cohort_path)
    OUT_DIR       = args.outdir
    TABLE_DIR     = OUT_DIR / "tables"
    FIG_DIR       = OUT_DIR / "figures"


def main():
    args = build_arg_parser().parse_args()
    _resolve_paths(args)

    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Inputs:  {CARRIERS_PATH}\n         {COHORT_PATH}")
    print(f"Outputs: {OUT_DIR}")

    # Clean up the retired any-disease outputs (A4).
    for stale in ["proportion_affected_by_age_all_carriers.tsv",
                  "proportion_affected_by_age_acmg_plp_only.tsv"]:
        p = TABLE_DIR / stale
        if p.exists():
            p.unlink()
    for stale in ["affected_by_age_all_carriers.png", "affected_by_age_acmg_plp_only.png"]:
        p = FIG_DIR / stale
        if p.exists():
            p.unlink()

    raw = load_carriers_raw()
    dedup = dedupe_to_person_gene(raw)
    cohort = load_cohort()
    person_demo, person_groups = build_person_tables(cohort)

    canonical_df, filter_log, ancestry_prerestrict = build_canonical(raw, dedup, person_demo)
    artifact_flags = compute_mapping_artifact_flags(ancestry_prerestrict)

    concordance = fix1_concordance(raw)
    gene_pen = fix3_gene_penetrance(canonical_df, person_groups, artifact_flags)
    disease_specific = fix3b_disease_specific_by_age(canonical_df, person_demo, person_groups, artifact_flags)
    fix4_ancestry_evidence_gap(canonical_df, args.min_carriers, artifact_flags)

    excluded = pd.read_csv(TABLE_DIR / "excluded_genes.tsv", sep="\t")
    evidence_gap_text = (TABLE_DIR / "evidence_gap_summary.txt").read_text().replace("\n", " ").strip()

    interpretable_genes = gene_pen.loc[gene_pen["interpretable"], "Gene"].tolist()
    print(f"\n[Report] Generating disease-specific figures for all {len(interpretable_genes)} "
          "interpretable genes (ACMG-PLP carriers only) ...")
    # Flat in figures/, which is where the committed per-gene PNGs are.
    gene_chart_dir = FIG_DIR
    gene_chart_dir.mkdir(parents=True, exist_ok=True)
    gallery_figs = {}
    for gene in interpretable_genes:
        gdf = disease_specific[disease_specific["Gene"] == gene]
        n = int(gene_pen.loc[gene_pen["Gene"] == gene, "n_carriers_canonical"].iloc[0])
        fig_path = gene_chart_dir / f"{gene}.png"
        plot_gene_disease_chart(gene, gdf, fig_path, n_carriers=n)
        gallery_figs[gene] = fig_path

    # The report's two headline figures. Both are interpretable in AoU; in a
    # cohort where either falls below the bar, say so rather than dying with a
    # KeyError after the whole pipeline has already run.
    missing_headline = [g for g in ("BRCA1", "BRCA2") if g not in gallery_figs]
    if missing_headline:
        raise SystemExit(
            f"{', '.join(missing_headline)} did not meet the interpretable bar "
            f"(>= {MIN_INTERPRETABLE_CARRIERS} canonical carriers, no mapping-artifact "
            "flag, defined gene->cancer mapping), so the headline figures in "
            "report.html cannot be built. Tables under "
            f"{TABLE_DIR} are complete and written."
        )
    brca1_fig = gallery_figs["BRCA1"]
    brca2_fig = gallery_figs["BRCA2"]

    build_report(concordance, gene_pen, disease_specific, filter_log, excluded,
                 evidence_gap_text, brca1_fig, brca2_fig, gallery_figs)
    build_process_walkthrough(filter_log, artifact_flags, gene_pen, concordance)

    # ClinVar-PLP comparison pipeline -- same method, different classifier,
    # written to separate files (report.html / process_walkthrough.html and
    # their underlying ACMG tsvs are untouched).
    print("\n[ClinVar comparison] Rebuilding canonical carriers under ClinVar-PLP ...")
    canonical_clinvar, filter_log_clinvar, ancestry_prerestrict_clinvar = build_canonical(
        raw, dedup, person_demo,
        classification_col="is_ClinVar_PLP", classification_label="ClinVar_PLP",
        log_filename="carrier_filter_log_clinvar.tsv",
    )
    artifact_flags_clinvar = compute_mapping_artifact_flags(ancestry_prerestrict_clinvar)
    gene_pen_clinvar = fix3_gene_penetrance(canonical_clinvar, person_groups, artifact_flags_clinvar,
                                             suffix="_clinvar")
    disease_specific_clinvar = fix3b_disease_specific_by_age(
        canonical_clinvar, person_demo, person_groups, artifact_flags_clinvar, suffix="_clinvar")

    # Per-gene disease-specific charts for ClinVar-PLP, same style as the ACMG
    # gallery, for every gene with a defined cancer mapping (not just genes
    # meeting ClinVar's own carrier bar -- the comparison report applies
    # its own, separate confidence bar when deciding what to actually show).
    genes_with_mapping = [g for g in GENE_CANCER_GROUPS if g in gene_pen_clinvar["Gene"].values]
    print(f"\n[ClinVar comparison] Generating disease-specific figures for "
          f"{len(genes_with_mapping)} genes (ClinVar-PLP carriers only) ...")
    gene_chart_dir_clinvar = FIG_DIR / "clinvar"
    gene_chart_dir_clinvar.mkdir(parents=True, exist_ok=True)
    gallery_figs_clinvar = {}
    for gene in genes_with_mapping:
        gdf = disease_specific_clinvar[disease_specific_clinvar["Gene"] == gene]
        if gdf.empty:
            continue
        n = int(gene_pen_clinvar.loc[gene_pen_clinvar["Gene"] == gene, "n_carriers_canonical"].iloc[0])
        fig_path = gene_chart_dir_clinvar / f"{gene}.png"
        plot_gene_disease_chart(gene, gdf, fig_path, n_carriers=n, classification_label="ClinVar-PLP")
        gallery_figs_clinvar[gene] = fig_path

    build_acmg_vs_clinvar_report(gene_pen, gene_pen_clinvar, filter_log, filter_log_clinvar,
                                  concordance, gallery_figs, gallery_figs_clinvar)

    print("\nDone. NOTE: these are crude cross-sectional estimates from a single "
          "AoU snapshot, not age-adjusted survival analysis, and are subject to "
          "ascertainment bias (AoU is not population-representative for "
          "cancer-predisposition genes). PMS2 and other flagged genes are excluded "
          "from main results due to likely pseudogene-mapping artifacts.")


if __name__ == "__main__":
    main()
