import pytest

ARRIBA_COLUMNS = [
    "#gene1", "gene2", "breakpoint1", "breakpoint2", "type",
    "split_reads1", "split_reads2", "discordant_mates", "confidence", "read_identifiers",
]
STARFUSION_COLUMNS = [
    "#FusionName", "JunctionReadCount", "SpanningFragCount", "LeftBreakpoint",
    "RightBreakpoint", "JunctionReads", "SpanningFrags", "FFPM",
]


def arriba_row(g1="TMPRSS2", g2="ERG", bp1="chr21:42870046", bp2="chr21:39817544",
               split1=10, split2=10, span=5, confidence="high", reads="."):
    return [g1, g2, bp1, bp2, "deletion", split1, split2, span, confidence, reads]


def starfusion_row(name="TMPRSS2--ERG", junction=10, span=5, left="chr21:42870046:-",
                   right="chr21:39817544:-", junction_reads=".", span_reads="."):
    return [name, junction, span, left, right, junction_reads, span_reads, "1.5"]


def write_tsv(path, columns, rows):
    lines = ["\t".join(columns)] + ["\t".join(str(v) for v in row) for row in rows]
    path.write_text("\n".join(lines) + "\n")


@pytest.fixture
def merge(run_script, read_tsv, tmp_path):
    def _merge(arriba_rows=(), starfusion_rows=None, extra=()):
        arriba = tmp_path / "a.tsv"
        write_tsv(arriba, ARRIBA_COLUMNS, arriba_rows)
        args = ["--sample", "S1", "--arriba", arriba,
                "--merged", tmp_path / "m.tsv", "--reported", tmp_path / "r.tsv", *extra]
        if starfusion_rows is not None:
            starfusion = tmp_path / "s.tsv"
            write_tsv(starfusion, STARFUSION_COLUMNS, starfusion_rows)
            args += ["--starfusion", starfusion]
        result = run_script("merge_fusions.py", *args)
        assert result.returncode == 0, result.stderr
        return read_tsv(tmp_path / "m.tsv"), read_tsv(tmp_path / "r.tsv")

    return _merge


def test_second_row_for_a_gene_pair_does_not_overwrite_the_first(merge):
    merged, reported = merge([
        arriba_row(bp1="chr21:100", split1=100, split2=100, confidence="high"),
        arriba_row(bp1="chr21:5000", split1=1, split2=0, span=0, confidence="low"),
    ])
    assert len(merged) == 2
    assert [r["left_breakpoint"] for r in reported] == ["chr21:100"]
    assert reported[0]["arriba_confidence"] == "high"


def test_gene_order_is_preserved(merge):
    merged, _ = merge([arriba_row()])
    assert merged[0]["fusion"] == "TMPRSS2--ERG"


def test_reciprocal_fusions_stay_separate(merge):
    merged, _ = merge([
        arriba_row(),
        arriba_row(g1="ERG", g2="TMPRSS2", bp1="chr21:39817544", bp2="chr21:42870046"),
    ])
    assert sorted(r["fusion"] for r in merged) == ["ERG--TMPRSS2", "TMPRSS2--ERG"]


def test_callers_agree_on_the_same_breakpoint(merge):
    merged, _ = merge([arriba_row()], [starfusion_row()])
    assert len(merged) == 1
    assert merged[0]["n_callers"] == "2"
    assert merged[0]["callers"] == "arriba,starfusion"


def test_strand_suffix_does_not_prevent_agreement(merge):
    merged, _ = merge([arriba_row(bp1="chr1:100", bp2="chr2:200")],
                      [starfusion_row(left="chr1:100:+", right="chr2:200:-")])
    assert len(merged) == 1


def test_breakpoints_within_window_agree(merge):
    merged, _ = merge([arriba_row(bp1="chr1:100", bp2="chr2:200")],
                      [starfusion_row(left="chr1:110:+", right="chr2:190:-")])
    assert len(merged) == 1
    assert merged[0]["n_callers"] == "2"


def test_breakpoints_outside_window_are_separate_events(merge):
    merged, _ = merge([arriba_row(bp1="chr1:100", bp2="chr2:200")],
                      [starfusion_row(left="chr1:111:+", right="chr2:200:-")])
    assert len(merged) == 2
    assert {r["n_callers"] for r in merged} == {"1"}


def test_window_is_configurable(merge):
    merged, _ = merge([arriba_row(bp1="chr1:100", bp2="chr2:200")],
                      [starfusion_row(left="chr1:111:+", right="chr2:200:-")],
                      extra=["--breakpoint-window", "11"])
    assert len(merged) == 1


def test_shared_reads_are_counted_once(merge):
    merged, _ = merge([arriba_row(split1=2, split2=1, span=0, reads="r1,r2,r3")],
                      [starfusion_row(junction=3, span=1, junction_reads="r1,r2,r3",
                                      span_reads="r4")])
    assert merged[0]["total_support"] == "4"
    assert merged[0]["support_source"] == "unique_read_names"
    assert merged[0]["arriba_split"] == "3"
    assert merged[0]["sf_split"] == "3"
    assert merged[0]["sf_span"] == "1"


def test_read_name_formatting_differences_do_not_inflate_support(merge):
    merged, _ = merge([arriba_row(split1=1, split2=0, span=0, reads="r1")],
                      [starfusion_row(junction=1, span=0, junction_reads="@r1/1")])
    assert merged[0]["total_support"] == "1"


def test_missing_read_names_fall_back_to_the_larger_caller(merge):
    merged, _ = merge([arriba_row(split1=10, split2=10, span=5, reads=".")],
                      [starfusion_row(junction=30, span=2)])
    assert merged[0]["total_support"] == "32"
    assert merged[0]["support_source"] == "caller_max"


def test_very_long_read_name_fields_do_not_crash(merge):
    names = ",".join(f"read_{i:06d}_{'x' * 40}" for i in range(6000))
    merged, _ = merge([arriba_row(split1=3000, split2=3000, span=0, reads=names)])
    assert merged[0]["total_support"] == "6000"
    assert merged[0]["n_read_names"] == "6000"


def test_read_names_in_output_are_capped(merge):
    names = ",".join(f"r{i}" for i in range(80))
    merged, _ = merge([arriba_row(split1=40, split2=40, span=0, reads=names)])
    assert len(merged[0]["read_names"].split(",")) == 50
    assert merged[0]["n_read_names"] == "80"


def test_empty_callers_give_valid_empty_outputs(merge):
    merged, reported = merge([], [])
    assert merged == []
    assert reported == []


def test_completely_empty_starfusion_file_is_valid(run_script, read_tsv, tmp_path):
    arriba = tmp_path / "a.tsv"
    write_tsv(arriba, ARRIBA_COLUMNS, [arriba_row()])
    starfusion = tmp_path / "s.tsv"
    starfusion.write_text("")
    result = run_script("merge_fusions.py", "--sample", "S1", "--arriba", arriba,
                        "--starfusion", starfusion, "--merged", tmp_path / "m.tsv",
                        "--reported", tmp_path / "r.tsv")
    assert result.returncode == 0, result.stderr
    assert len(read_tsv(tmp_path / "m.tsv")) == 1


def test_intergenic_gene_names_are_kept_verbatim(merge):
    merged, _ = merge([arriba_row(g1="A(100),B(200)", g2="C", bp1="chr1:1", bp2="chr2:2")])
    assert merged[0]["fusion"] == "A(100),B(200)--C"


def test_two_callers_promote_a_low_confidence_call(merge):
    _, reported = merge([arriba_row(confidence="low")], [starfusion_row()])
    assert len(reported) == 1


def test_single_caller_medium_confidence_is_not_promoted(merge):
    _, reported = merge([arriba_row(confidence="medium")])
    assert reported == []


def test_rows_are_ranked_by_callers_then_support(merge):
    merged, _ = merge(
        [
            arriba_row(g1="A", g2="B", bp1="chr1:1", bp2="chr2:1", split1=90, split2=90),
            arriba_row(g1="C", g2="D", bp1="chr3:1", bp2="chr4:1", split1=1, split2=1, span=0),
        ],
        [starfusion_row(name="C--D", junction=2, span=0, left="chr3:1:+", right="chr4:1:+")],
    )
    assert [r["fusion"] for r in merged] == ["C--D", "A--B"]
