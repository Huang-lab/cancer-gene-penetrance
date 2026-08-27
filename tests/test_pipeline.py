"""End-to-end checks on the pipeline, run against a synthetic cohort.

`test_outputs_land_under_outdir` is the regression test for this PR: the script
used to write its tables and HTML into the repo root and its figures into
`./figures/`, while reading its inputs from one level *above* the repo. Nothing
under `results/` was reachable by running it, so README.md's "Running the script
regenerates everything under results/" was not true of any invocation.
"""
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent))
import synthetic_cohort  # noqa: E402

REPO = Path(__file__).parent.parent
SCRIPT = REPO / "gene_penetrance.py"

# The exact file names the committed results/ tree contains.
EXPECTED_TABLES = {
    "acmg_clinvar_concordance_by_ancestry.tsv",
    "carrier_filter_log.tsv",
    "carrier_filter_log_clinvar.tsv",
    "disease_specific_penetrance_by_age.tsv",
    "disease_specific_penetrance_by_age_clinvar.tsv",
    "evidence_gap_summary.txt",
    "excluded_genes.tsv",
    "excluded_genes_clinvar.tsv",
    "gene_by_ancestry_counts.tsv",
    "gene_penetrance.tsv",
    "gene_penetrance_clinvar.tsv",
}
EXPECTED_HTML = {"report.html", "process_walkthrough.html", "report_acmg_vs_clinvar.html"}


@pytest.fixture(scope="session")
def run(tmp_path_factory):
    """Run the whole pipeline once on synthetic data; reuse for every test."""
    indir = tmp_path_factory.mktemp("in")
    outdir = tmp_path_factory.mktemp("out") / "results"
    carriers, cohort = synthetic_cohort.write(indir)

    proc = subprocess.run(
        [sys.executable, str(SCRIPT),
         "--carriers-path", str(carriers),
         "--cohort-path", str(cohort),
         "--outdir", str(outdir)],
        capture_output=True, text=True, cwd=REPO,
    )
    assert proc.returncode == 0, f"pipeline failed:\n{proc.stdout}\n{proc.stderr}"
    return outdir


def test_outputs_land_under_outdir(run):
    """Tables in tables/, figures in figures/, reports at the top -- the layout
    that is committed. Before this PR every one of these went to the repo root."""
    assert EXPECTED_TABLES <= {p.name for p in (run / "tables").iterdir()}
    assert EXPECTED_HTML <= {p.name for p in run.glob("*.html")}
    assert (run / "figures" / "concordance_by_ancestry.png").exists()
    assert (run / "figures" / "BRCA1.png").exists(), "per-gene charts sit flat in figures/"

    # Nothing may be written beside the script any more.
    strays = EXPECTED_TABLES | EXPECTED_HTML
    assert not [n for n in strays if (REPO / n).exists()], "output written to the repo root"


def test_missing_input_is_reported_clearly(tmp_path):
    """A wrong --carriers-path should name the file, not surface as a pandas
    traceback from inside load_carriers_raw."""
    proc = subprocess.run(
        [sys.executable, str(SCRIPT),
         "--carriers-path", str(tmp_path / "nope.tsv"),
         "--cohort-path", str(tmp_path / "also-nope.csv.gz"),
         "--outdir", str(tmp_path / "out")],
        capture_output=True, text=True, cwd=REPO,
    )
    assert proc.returncode != 0
    assert "no such file" in proc.stderr.lower()
    assert "nope.tsv" in proc.stderr
    assert "Traceback" not in proc.stderr


def test_gene_below_the_bar_is_excluded_with_a_reason(run):
    """The minimum-evidence bar is documented as never dropping a gene silently.
    WT1 is seeded below it in the synthetic cohort."""
    excluded = pd.read_csv(run / "tables" / "excluded_genes.tsv", sep="\t")
    wt1 = excluded[excluded["Gene"] == "WT1"]
    assert len(wt1) == 1
    assert "canonical carriers" in wt1.iloc[0]["exclusion_reason"]

    penetrance = pd.read_csv(run / "tables" / "gene_penetrance.tsv", sep="\t")
    assert not penetrance.loc[penetrance["Gene"] == "WT1", "interpretable"].iloc[0]


def test_carrier_filter_log_accounts_for_every_drop(run):
    """The canonical denominator's selling point is that each step is counted.
    Each stage must be monotonically non-increasing, and the cohort-join drop
    must equal the difference it claims to explain."""
    log = pd.read_csv(run / "tables" / "carrier_filter_log.tsv", sep="\t")
    for _, r in log.iterrows():
        assert r["unique_person_gene_ACMG_PLP"] <= r["unique_person_gene_any_classification"]
        assert r["joined_to_cohort"] <= r["unique_person_gene_ACMG_PLP"]
        assert r["n_carriers_canonical"] <= r["joined_to_cohort"]
        assert (r["dropped_at_cohort_join"]
                == r["unique_person_gene_ACMG_PLP"] - r["joined_to_cohort"])


@pytest.mark.xfail(
    strict=True,
    reason="Sex restriction is applied per gene, not per (gene, cancer). See "
           "issue: 'The sex restriction is applied per gene, so Lynch colorectal "
           "penetrance is computed on women only'.",
)
def test_one_disease_has_one_background_rate(run):
    """A population background rate is a property of (disease, age band) and the
    cohort -- not of which gene's row you happen to be reading.

    It currently is not. `FEMALE_ONLY_GENES` restricts the whole gene, so the
    Colorectum baseline is computed over women only for MLH1/MSH2 and over
    everyone for APC, and both land in the same column of the same table.
    """
    ds = pd.read_csv(run / "tables" / "disease_specific_penetrance_by_age.tsv", sep="\t")
    for (disease, band), grp in ds.groupby(["disease", "age_band"]):
        rates = grp["baseline_rate"].dropna().unique()
        assert len(rates) <= 1, (
            f"{disease} / {band} has {len(rates)} different background rates "
            f"({rates}) depending on the gene: "
            f"{sorted(grp['Gene'].unique())}"
        )
