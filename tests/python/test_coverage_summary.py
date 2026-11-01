import gzip
import shutil
import subprocess

import pytest

FLAGSTAT = """\
1200 + 0 in total (QC-passed reads + QC-failed reads)
1000 + 0 primary
150 + 0 secondary
50 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
1100 + 0 mapped (91.67% : N/A)
900 + 0 primary mapped (90.00% : N/A)
1000 + 0 paired in sequencing
"""


@pytest.fixture
def summarise(run_script, read_tsv, tmp_path):
    def _summarise(sample="S1", flagstat=FLAGSTAT):
        thresholds = tmp_path / f"{sample}.thresholds.bed.gz"
        with gzip.open(thresholds, "wt") as fh:
            fh.write("#chrom\tstart\tend\tregion\t50X\t100X\t200X\t500X\t1000X\n")
            fh.write("chr1\t0\t100\tA\t100\t50\t0\t0\t0\n")
        summary = tmp_path / f"{sample}.summary.txt"
        summary.write_text("chrom\tlength\tbases\tmean\tmin\tmax\n"
                           "total_region\t100\t5000\t50.00\t0\t90\n")
        stats = tmp_path / f"{sample}.flagstat"
        stats.write_text(flagstat)
        out = tmp_path / f"{sample}.tsv"
        mqc = tmp_path / f"{sample}.coverage_mqc.txt"
        result = run_script("coverage_summary.py", "--sample", sample, "--thresholds", thresholds,
                            "--summary", summary, "--flagstat", stats, "--out", out, "--mqc", mqc)
        assert result.returncode == 0, result.stderr
        return read_tsv(out)[0], mqc

    return _summarise


def test_total_reads_excludes_secondary_and_supplementary(summarise):
    row, _ = summarise()
    assert row["total_reads"] == "1000"
    assert row["total_read_pairs"] == "500"


def test_old_flagstat_without_primary_line_is_derived(summarise):
    row, _ = summarise(flagstat="120 + 0 in total (QC-passed reads + QC-failed reads)\n"
                                "15 + 0 secondary\n5 + 0 supplementary\n")
    assert row["total_reads"] == "100"
    assert row["total_read_pairs"] == "NA"


def test_percentages_are_calculated(summarise):
    row, _ = summarise()
    assert row["pct_bases_ge_50x"] == "100.0"
    assert row["pct_bases_ge_100x"] == "50.0"


def test_multiqc_shows_real_samples_and_no_phantom_samples(summarise, tmp_path):
    if not shutil.which("multiqc"):
        pytest.skip("multiqc is not installed")
    summarise("S1")
    summarise("S2")
    outdir = tmp_path / "report"
    subprocess.run(["multiqc", "-q", "-f", "-o", str(outdir), str(tmp_path)], check=True)
    general_stats = (outdir / "multiqc_data" / "multiqc_general_stats.txt").read_text()
    samples = {line.split("\t")[0] for line in general_stats.splitlines()[1:]}
    assert samples == {"S1", "S2"}
    assert "total_reads" not in samples
    assert "mean_target_coverage" in general_stats.splitlines()[0]
