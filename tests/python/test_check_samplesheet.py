import pytest

HEADER = "sample,fastq_1,fastq_2,strandedness\n"


def check(run_script, tmp_path, content):
    infile = tmp_path / "in.csv"
    outfile = tmp_path / "out.csv"
    infile.write_bytes(content.encode("utf-8"))
    return run_script("check_samplesheet.py", infile, outfile), outfile


def test_extra_field_gives_clear_error_naming_the_line(run_script, tmp_path):
    result, _ = check(run_script, tmp_path, HEADER + "S1,a.fastq.gz,b.fastq.gz,reverse,extra\n")
    assert result.returncode == 1
    assert "line 2" in result.stderr
    assert "Traceback" not in result.stderr


def test_byte_order_mark_is_accepted(run_script, tmp_path):
    result, outfile = check(
        run_script, tmp_path, "﻿" + HEADER + "S1,a.fastq.gz,b.fastq.gz,reverse\n"
    )
    assert result.returncode == 0, result.stderr
    assert outfile.read_text().splitlines()[1].startswith("S1,")


def test_identical_read_files_are_rejected(run_script, tmp_path):
    result, _ = check(run_script, tmp_path, HEADER + "S1,a.fastq.gz,a.fastq.gz,reverse\n")
    assert result.returncode == 1
    assert "line 2" in result.stderr


@pytest.mark.parametrize("sample", ["S/1", "S 1", "S;1", "../S1", "S1*"])
def test_unsafe_sample_ids_are_rejected(run_script, tmp_path, sample):
    result, _ = check(run_script, tmp_path, HEADER + f"{sample},a.fastq.gz,b.fastq.gz,reverse\n")
    assert result.returncode == 1
    assert "line 2" in result.stderr


@pytest.mark.parametrize(
    "row, expected",
    [
        ("S1,a.fastq.gz,b.fastq.gz,reverse", "S1,a.fastq.gz,b.fastq.gz,reverse"),
        ("S.1-a_b,a.fq.gz,,forward", "S.1-a_b,a.fq.gz,,forward"),
        ("S2,a.fq.gz,,", "S2,a.fq.gz,,reverse"),
    ],
)
def test_valid_rows_pass_through(run_script, tmp_path, row, expected):
    result, outfile = check(run_script, tmp_path, HEADER + row + "\n")
    assert result.returncode == 0, result.stderr
    assert outfile.read_text().splitlines() == [HEADER.strip(), expected]
