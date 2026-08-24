# Transcript — ACMG vs ClinVar Comparison (report_acmg_vs_clinvar.html + process_walkthrough.html)

## The whole picture

The main report defines "carrier" using ACMG-PLP classification. That choice matters enough to
justify checking directly: what would the same analysis have produced if we'd required ClinVar-PLP
confirmation instead — the more common default in a lot of clinical pipelines? This report answers
that by rerunning the identical pipeline a second time, changing exactly one input, so the
difference between the two outputs is attributable only to the classifier choice.

## What we did

1. **Reused every step of the main pipeline unchanged** — same 27 genes, same gene→cancer-type
   mapping, same canonical-carrier-denominator method, same background-rate comparison — and swapped
   only the restriction step: `is_ACMG_PLP == True` became `is_ClinVar_PLP == True`. Everything else,
   including the sex restriction and the cohort join, stayed identical.

2. **Compared carrier counts and penetrance gene-by-gene.** For every gene with a defined cancer
   mapping, we now have two canonical carrier counts (ACMG N, ClinVar N) and two penetrance estimates
   side by side.

3. **Applied one confidence bar to both sides.** A gene only gets shown as a figure if ACMG N ≥ 20
   AND ClinVar N ≥ 20 AND it isn't flagged as a pseudogene artifact — the same threshold used on both
   classifiers, so neither side gets an easier standard. Genes that fail are listed with the specific
   reason (which side, and by how much) rather than dropped from the underlying tables.

4. **Generated matched per-gene charts for both classifiers** — same chart style as the main report
   (each cancer type its own line, population background dashed, ★ for significant elevation) — so
   the ACMG and ClinVar curves for the same gene can be viewed side by side directly, not just as
   summary numbers.

## Major results

**Aggregate cost.** Across the 27 genes: ACMG-PLP carrier total = 3,790. ClinVar-PLP carrier total =
2,066. Requiring ClinVar confirmation drops the usable carrier pool by **45.5%** — this is the
concrete, carrier-count version of the abstract discordance percentages in the main report's
Finding 1.

**Which genes lose the most.** TSC2 keeps only 2.3% of its ACMG carriers under ClinVar (306 → 7).
SMAD4 keeps 2.2% (134 → 3). MAX keeps 5.1% (39 → 2). VHL keeps 14.4% (104 → 15). These genes fall
below the 20-carrier bar on the ClinVar side and are excluded from the figure gallery as a direct
consequence — their absence from the chart IS the finding.

**Which genes lose the least.** SDHC retains 84.5% (58 → 49), SDHB 82.2%, BRCA2 80.0% (705 → 564).

**BRCA1 is a genuine exception, worth naming explicitly.** ACMG gives 320 carriers (44.7%
penetrance); ClinVar gives 381 carriers (46.7% penetrance) — ClinVar actually classifies *more*
BRCA1 variants as pathogenic overall than our ACMG ruleset does. This doesn't undercut the ancestry
argument — the concordance finding is about which carriers get dropped disproportionately by
ancestry, not about every single gene's raw total moving in the same direction. BRCA2 shows the more
typical pattern: ACMG 705 (29.4%) vs. ClinVar 564 (27.8%) — smaller gap, same direction.

**Figure coverage.** 15 of 27 genes clear the confidence bar and get paired ACMG-vs-ClinVar charts.
12 are excluded: TSC2, SMAD4, MAX, SDHAF2, TSC1, STK11, RB1, NF2, WT1, BMPR1A (ClinVar N too low),
VHL (ClinVar N=15, just under 20), and PMS2 (pseudogene artifact, excluded on both sides
consistently with the main report).

## Statistics used

Same as the main report: 95% Wilson score interval on every proportion (penetrance, carrier-count
ratio is a plain ratio, not a CI). "Significantly elevated" = carrier-rate CI lower bound > background
CI upper bound, non-overlapping-interval test, no p-value computed. Confidence bar for figures (ACMG
N≥20 AND ClinVar N≥20, no artifact flag) is a fixed inclusion rule, not a statistical test. Full
formulas in report_acmg_vs_clinvar.html under "Calculation and statistics used."

## Interpretation

This is not a "which classifier is right" comparison — ACMG-PLP is the definition used throughout
the main report, and the reasoning stands: it's ancestry-neutral, ClinVar measurably isn't. What this
report adds is a concrete number for what that choice costs if you'd gone the other way — a 45.5%
carrier loss, concentrated in specific genes, in the direction the ancestry-discordance numbers
already predicted. It turns Finding 1 from a percentage into a carrier count anyone can point to.
