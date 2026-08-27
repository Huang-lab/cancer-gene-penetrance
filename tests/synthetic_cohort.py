"""Generate a synthetic carrier + cohort pair matching README.md's input schema.

No participant-level data ships with this repo, so every test here runs against
data made up on the spot. It is shaped to exercise the pipeline's branches --
both sexes, several ancestries, sex-specific and non-sex-specific cancers, genes
above and below the interpretable bar -- not to look like real epidemiology.
"""
import gzip
from pathlib import Path

import numpy as np
import pandas as pd

ANCESTRIES = ["eur", "afr", "amr", "sas", "eas"]

# (gene, n_carriers, fraction of carriers given the gene's first mapped cancer)
GENE_SPEC = [
    ("BRCA1", 60, 0.40),
    ("BRCA2", 60, 0.30),
    ("MLH1", 60, 0.25),
    ("MSH2", 40, 0.15),
    ("APC", 40, 0.15),
    ("TP53", 30, 0.10),
    ("PTEN", 30, 0.10),
    ("SDHB", 25, 0.04),
    ("WT1", 8, 0.10),   # below the interpretable bar on purpose
]

GENE_FIRST_CANCER = {
    "BRCA1": "Breast", "BRCA2": "Breast", "MLH1": "Colorectum", "MSH2": "Colorectum",
    "APC": "Colorectum", "TP53": "Breast", "PTEN": "Breast", "SDHB": "Endocrine System",
    "WT1": "Kidney",
}


def write(tmpdir: Path, seed: int = 0) -> tuple:
    """Write carriers.tsv and cohort.csv.gz into tmpdir; return their paths."""
    rng = np.random.default_rng(seed)
    carrier_rows, cohort_rows = [], []
    pid = 0

    # A background population so the per-disease baseline rates have a denominator.
    for _ in range(4000):
        pid += 1
        cohort_rows.append({
            "person_id": pid,
            "age": float(rng.integers(30, 85)),
            "sex_at_birth": rng.choice(["Female", "Male"]),
            "Group": rng.choice(
                ["Control (age >= 50)", "Control (age < 50)", "Breast", "Colorectum",
                 "Ovary", "Endometrium", "Kidney", "Endocrine System", "Thyroid"],
                p=[0.40, 0.30, 0.05, 0.05, 0.04, 0.04, 0.04, 0.04, 0.04],
            ),
        })

    for gene, n, affected_frac in GENE_SPEC:
        for i in range(n):
            pid += 1
            sex = "Female" if i % 2 == 0 else "Male"
            age = float(rng.integers(30, 85))
            affected = rng.random() < affected_frac

            # Most carriers get one variant; some get two in the same gene, which is
            # what makes the person x gene collapse rule observable.
            n_variants = 2 if i % 7 == 0 else 1
            for v in range(n_variants):
                carrier_rows.append({
                    "person_id": pid,
                    "Gene": gene,
                    "Variant_ID": f"chr1:{100000 + i * 10 + v}:A>G",
                    "is_ACMG_PLP": True,
                    # Deliberately mixed within a person: variant 0 confirmed in
                    # ClinVar, variant 1 not.
                    "is_ClinVar_PLP": (v == 0),
                    "ancestry_pred": ANCESTRIES[i % len(ANCESTRIES)],
                    "sex_at_birth": sex,
                    "age": age,
                })

            cohort_rows.append({
                "person_id": pid, "age": age, "sex_at_birth": sex,
                "Group": GENE_FIRST_CANCER[gene] if affected else "Control (age >= 50)",
            })

    carriers = pd.DataFrame(carrier_rows)
    cohort = pd.DataFrame(cohort_rows)

    carriers_path = tmpdir / "carriers.tsv"
    cohort_path = tmpdir / "cohort.csv.gz"
    carriers.to_csv(carriers_path, sep="\t", index=False)
    with gzip.open(cohort_path, "wt") as fh:
        cohort.to_csv(fh, index=False)
    return carriers_path, cohort_path
