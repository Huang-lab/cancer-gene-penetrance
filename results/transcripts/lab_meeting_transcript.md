# Transcript — Main Results (report.html + process_walkthrough.html)

## The whole picture

We have AoU carrier data for 27 cancer-predisposition genes (BRCA1, BRCA2, Lynch-syndrome genes,
paraganglioma genes, etc.) and cohort data linking each person to age, sex, ancestry, and any
recorded cancer diagnosis. Two things are worth measuring from that: whether carrying a
likely-pathogenic variant actually predicts getting the matched cancer (penetrance), and whether the
tools we use to call a variant "likely-pathogenic" work the same way regardless of a person's
ancestry (concordance). Both come out of the same pipeline, run once, applied consistently across
all 27 genes.

## What we did

1. **Built one canonical carrier list per gene.** Earlier versions of this analysis had three
   different carrier counts for BRCA1 in three different tables (565, 611, 917) because each table
   quietly built its own denominator. We fixed this by defining a single fixed sequence — start from
   every raw variant row for the gene, collapse to one entry per person, keep only people whose
   variant is called ACMG-PLP (pathogenic/likely-pathogenic by rule-based classification), join to
   the cohort table to get their phenotype (dropping anyone with no matching cohort record, and
   counting that drop explicitly), then restrict to the correct sex if the gene's cancer type is
   sex-specific (e.g. BRCA1 restricted to females). Whatever number survives that chain is used,
   unchanged, everywhere that gene appears.

2. **Chose ACMG-PLP, not ClinVar-PLP, as the carrier definition.** ACMG is a rule-based classifier
   applied the same way no matter who the patient is. ClinVar is a database of manually submitted
   classifications that is measurably less complete for non-European variants (that's Finding 1
   below) — so defining "carrier" by ClinVar confirmation would have baked that same ancestry bias
   into every penetrance number before we even started counting.

3. **Computed penetrance per gene, per age band, against a real background.** For each gene we use an
   explicit gene→cancer-type map (BRCA1→breast/ovary, MLH1→colorectal/endometrial, etc.) rather than
   counting any cancer as a "hit" — a BRCA1 carrier with lung cancer isn't a BRCA1 penetrance event.
   For every (gene, cancer type, age band) we compute: carriers with that cancer ÷ all carriers in
   that age band, with a 95% confidence interval. We then compute the same disease's rate in the
   *entire* cohort (not just carriers) in that same age band and sex, as a background comparison —
   and only mark a carrier rate "significantly elevated" when its confidence interval doesn't overlap
   the background's.

4. **Computed concordance per ancestry.** Among people whose variant is ACMG-PLP, what fraction are
   *not* also confirmed ClinVar-PLP — broken out by predicted ancestry.

5. **Flagged and excluded known artifacts.** PMS2 has a well-documented pseudogene mismapping problem
   (PMS2/PMS2CL) that inflates carrier counts in specific ancestries for technical, not biological,
   reasons — it's excluded from all main results. Other genes get the same exclusion if one non-EUR
   ancestry shows more than double the EUR carrier count with no plausible biological reason.

6. **Set a minimum-evidence bar.** A gene needs ≥20 canonical carriers, no artifact flag, and a
   defined cancer mapping to appear in the main results ("interpretable"). 23 of 27 genes clear that
   bar; the other 4 are listed separately with the specific reason, not silently dropped.

## Major results

**Finding 1 — concordance gap by ancestry.** Among ACMG-PLP carriers, the fraction with NO matching
ClinVar-PLP confirmation:

| Ancestry | Discordance |
|---|---|
| EUR | 42.5% |
| AFR | 51.9% |
| AMR | 53.2% |
| SAS | 67.8% |
| EAS | 90.9% |

Every non-EUR confidence interval sits clearly above EUR's — not noise. Practical meaning: a pipeline
that trusts ClinVar as the pathogenicity gate will under-ascertain real carriers far more heavily in
non-European populations than in European ones. One caveat we checked and are reporting explicitly:
the EAS 90.9% figure is 83.5% attributable to a single recurrent MSH6 variant
(chr2:47806842:TTTGA>TTTGATTGA) — it should be read as that one variant's classification status, not
as evidence of broad-based EAS under-curation. The other ancestries' gaps are spread across many
different variants, which is the more typical under-curation signature.

**Finding 2 — penetrance behaves as expected.**
- BRCA1: 320 canonical carriers, 44.7% develop breast or ovarian cancer, rising through the 40s/50s,
  above background at every age band.
- BRCA2: 705 carriers, 29.4% penetrance, same shape, consistently below BRCA1 — matches the known
  BRCA1 > BRCA2 risk ordering.
- All 23 interpretable genes have their own age-stratified chart in the report gallery, not just
  these two — BRCA1/BRCA2 lead because they have the largest N and tightest confidence intervals.
- Excluded from main results: BMPR1A (n=13), NF2 (n=19), WT1 (n=16) — below the 20-carrier bar;
  PMS2 (n=263) — pseudogene artifact.
- These numbers run lower than clinic-ascertained penetrance estimates you may have seen elsewhere
  (which select for high-risk families already flagged for testing) — expected for an unselected
  volunteer cohort, not a red flag.

## Statistics used

All proportions (penetrance, background rate, discordance rate) use a 95% **Wilson score interval**,
not the normal approximation — several cells have small N or proportions near 0/1, where the normal
approximation undercoverages and can produce out-of-range bounds. Computed via
`statsmodels.stats.proportion.proportion_confint(method="wilson")`. "Significantly elevated" (★) is a
non-overlapping-interval test — carrier rate's CI lower bound > background rate's CI upper bound —
not a formal two-proportion hypothesis test; no p-value is computed. Full formulas are in report.html
under "How the numbers are calculated."

## Limitations, stated plainly

Crude, cross-sectional proportions — not survival/cumulative-incidence estimates, no adjustment for
follow-up time. AoU is a healthcare-utilizing volunteer cohort, not a random population sample.
Roughly 30% of ACMG-PLP carriers per gene have no matching cohort/phenotype record and are dropped
before any calculation — reported explicitly in the filter log, not hidden. The pseudogene-artifact
rule is a heuristic (PMS2 always flagged, a 2×-EUR-ratio check for two watch-list genes), not an
exhaustive QC pass. Single-cohort, preliminary, not yet validated against an independent dataset.
