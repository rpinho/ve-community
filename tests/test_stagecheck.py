"""Invented metadata only; no challenge files, measurements or network access."""
import argparse
import csv
import io
from decimal import Decimal as D

import pytest

from ve_community.stagecheck import (
    Age,
    check_rows,
    choose_column,
    main,
    margin,
    parse_age,
    read_sheet,
    render,
)


def check(stage, task="t1", **kwargs):
    row = {"sample": "invented-a", "stage": stage}
    row.update(kwargs.pop("annotations", {}))
    return check_rows([row], task=task, sample_column="sample", stage_column="stage", **kwargs)[0]


@pytest.mark.parametrize("task,stage,expected", [
    ("t1", "E9.5", "OUTSIDE_BY_STAGE"),
    ("t1", "E9.500001", "EXCLUDE"),
    ("t1", "E10.5", "EXCLUDE"),
    ("t1", "E12.5", "EXCLUDE"),
    ("t1", "E13.5", "EXCLUDE"),
    ("t1", "E13.500001", "OUTSIDE_BY_STAGE"),
    ("t2-heart", "E8.25", "OUTSIDE_BY_STAGE"),
    ("t2-heart", "E8.250001", "EXCLUDE"),
    ("t2-heart", "E8.5", "EXCLUDE"),
    ("t2-heart", "E8.75", "OUTSIDE_BY_STAGE"),
    ("t2-heart", "E9.5", "OUTSIDE_BY_STAGE"),
    ("t2-heart", "E10.0", "EXCLUDE"),
    ("t2-heart", "E13.5", "EXCLUDE"),
    ("t2-heart", "E13.6", "OUTSIDE_BY_STAGE"),
    ("t2-embryo", "E7.25", "OUTSIDE_BY_STAGE"),
    ("t2-embryo", "E7.250001", "EXCLUDE"),
    ("t2-embryo", "E7.5", "EXCLUDE"),
    ("t2-embryo", "E7.75", "EXCLUDE"),
    ("t2-embryo", "E8.0", "OUTSIDE_BY_STAGE"),
])
def test_published_window_boundaries(task, stage, expected):
    assert check(stage, task).decision == expected


@pytest.mark.parametrize("text,lo,hi", [
    ("E7.25", "7.25", "7.25"), ("7.25 days", "7.25", "7.25"),
    ("7.25 dpc", "7.25", "7.25"), ("E7.25-E8.0", "7.25", "8.0"),
    ("E7.25 to E8.0", "7.25", "8.0"), ("7.25–8", "7.25", "8"),
    ("E8 ± 0.25", "7.75", "8.25"), ("E8 +/- 0.25 days", "7.75", "8.25"),
])
def test_explicit_age_grammar(text, lo, hi):
    assert parse_age(text) == Age(D(lo), D(hi))


@pytest.mark.parametrize("text", [
    "", "nan", "inf", "-1", "31", "P0", "TS12", "12 somites", "litter E8",
    "embryo-7.5", "E8/E9", "E8-7", "~E8.0", "E8.0 (approximate)", "E8.0 hours",
])
def test_ambiguous_age_never_inferred(text):
    assert parse_age(text) is None
    assert check(text).decision == "REVIEW"


@pytest.mark.parametrize("task,stage,relation,decision", [
    ("t1", "E9-E9.5", "OUTSIDE", "OUTSIDE_BY_STAGE"),
    ("t1", "E9.5-E13.5", "CROSSES", "REVIEW"),
    ("t1", "E10-E13.5", "INSIDE", "EXCLUDE"),
    ("t1", "E13.5-E14", "CROSSES", "REVIEW"),
    ("t1", "E13.6-E14", "OUTSIDE", "OUTSIDE_BY_STAGE"),
    ("t2-heart", "E8.25-E8.75", "CROSSES", "REVIEW"),
    ("t2-heart", "E8.75-E9.5", "OUTSIDE", "OUTSIDE_BY_STAGE"),
    ("t2-heart", "E8-E14", "CROSSES", "REVIEW"),
    ("t2-embryo", "E7.25-E8", "CROSSES", "REVIEW"),
    ("t2-embryo", "E7.3-E7.9", "INSIDE", "EXCLUDE"),
])
def test_age_interval_intersection_includes_interior(task, stage, relation, decision):
    verdict = check(stage, task)
    assert (verdict.nominal_relation, verdict.decision) == (relation, decision)


@pytest.mark.parametrize("annotations", [
    {"somite_count": "0"}, {"Theiler stage": "TS12"},
    {"age_basis": "litter"}, {"staging method": "litter date"},
    {"litter_age": "E9.5"}, {"genotype": "Gata4 KO"},
    {"Characteristics[organism]": "Homo sapiens"},
    {"Characteristics[developmental stage]": "12 somites"},
    {"Characteristics[developmental stage]": "TS12"},
    {"Characteristics[age unit]": "hour"},
    {"Characteristics[developmental stage]": "litter E9.5"},
])
def test_nearness_and_condition_annotations_force_review(annotations):
    assert check("E9.5", annotations=annotations).decision == "REVIEW"


def test_litter_identifier_is_not_an_age_assignment():
    assert check("E9.5", annotations={"litter_id": "invented-litter"}).decision == "OUTSIDE_BY_STAGE"


def test_missing_annotation_is_not_evidence_and_wild_type_is_explicit():
    assert check("E9.5", annotations={"somite_count": "NA", "genotype": "WT"}).decision == "OUTSIDE_BY_STAGE"


def test_alternative_age_columns_cannot_silently_disagree():
    assert check("E9.5", annotations={"age": "E10"}).decision == "REVIEW"
    assert check("E9.5", annotations={"age": "9.5 days"}).decision == "OUTSIDE_BY_STAGE"


def test_uncertainty_crosses_boundary_and_nearness_band_is_separate():
    assert check("E9.5", uncertainty=D("0.1")).nominal_relation == "CROSSES"
    verdict = check("E9.4", near=D("0.1"))
    assert (verdict.decision, verdict.nominal_relation) == ("REVIEW", "OUTSIDE")
    assert "not a rule" in " ".join(verdict.reasons)


@pytest.mark.parametrize("value", ["NaN", "Infinity", "-0.1", "bad"])
def test_bad_margin_rejected(value):
    with pytest.raises(argparse.ArgumentTypeError):
        margin(value)


def test_t3_never_approved_by_age_alone():
    for stage in ("E8.75", "E9.5", "E14"):
        assert check(stage, "t3").decision == "REVIEW"


def test_grouping_repeated_sdrf_rows_and_conflicting_specimens():
    rows = [
        {"sample": "invented-a", "stage": "E9"},
        {"sample": "invented-a", "stage": "E9.0"},
        {"sample": "invented-b", "stage": "E9"},
        {"sample": "invented-b", "stage": "E10"},
        {"sample": "invented-c", "stage": "E8"},
        {"sample": "invented-c", "stage": "E9"},
        {"sample": "invented-d", "stage": "E9"},
        {"sample": "invented-d", "stage": "unknown"},
    ]
    a, b, c, d = check_rows(rows, task="t1", sample_column="sample", stage_column="stage")
    assert a.rows == (1, 2) and a.decision == "OUTSIDE_BY_STAGE"
    assert b.rows == (3, 4) and b.decision == "EXCLUDE"
    assert c.decision == d.decision == "REVIEW"


def test_explicit_columns_and_ambiguous_headers():
    headers = ["Source Name", "Characteristics[age]", "Characteristics[developmental stage]"]
    assert choose_column(headers, "sample", None) == "Source Name"
    with pytest.raises(ValueError, match="ambiguous"):
        choose_column(headers, "stage", None)
    assert choose_column(headers, "stage", headers[1]) == headers[1]
    with pytest.raises(ValueError, match="absent"):
        choose_column(headers, "stage", "not-present")


@pytest.mark.parametrize("content", [
    "sample,stage,stage\na,E9,E9\n", "sample,stage\na\n", "sample,stage\na,E9,extra\n",
    "sample,stage\n", "sample,\na,E9\n", "sample,Sample,stage\na,b,E9\n",
])
def test_malformed_sheets_rejected(tmp_path, content):
    sheet = tmp_path / "invented.csv"
    sheet.write_text(content)
    with pytest.raises(ValueError):
        read_sheet(sheet)


def test_sdrf_quoted_fields_bom_and_duplicate_records(tmp_path):
    sheet = tmp_path / "invented.sdrf.txt"
    sheet.write_text('\ufeffSource Name\tCharacteristics[developmental stage]\tComment[run]\n'
                     '"invented\tname"\tE9\tinvented-run-a\n'
                     '"invented\tname"\tE9\tinvented-run-b\n')
    headers, rows = read_sheet(sheet)
    assert len(headers) == 3
    assert len(rows) == 2
    verdicts = check_rows(rows, task="t1", sample_column=headers[0], stage_column=headers[1])
    assert verdicts[0].rows == (1, 2)
    parsed = list(csv.reader(io.StringIO(render(verdicts, "t1")), delimiter="\t"))
    assert parsed[1][0] == "invented\tname"


def test_formula_like_identifiers_are_escaped():
    v = check_rows([{"sample": "=synthetic", "stage": "E9"}], task="t1",
                   sample_column="sample", stage_column="stage")
    assert render(v, "t1").splitlines()[1].startswith("'=synthetic\t")


def test_cli_end_to_end_with_synthetic_fixture(tmp_path, capsys):
    sheet = tmp_path / "invented.csv"
    sheet.write_text("sample,stage\ninvented-a,E9.5\ninvented-b,E13.5\n")
    before = sheet.read_bytes()
    assert main([str(sheet), "--task", "t1"]) == 1
    output = capsys.readouterr().out
    assert "invented-a" in output and "OUTSIDE_BY_STAGE" in output
    assert "invented-b" in output and "EXCLUDE" in output
    assert sheet.read_bytes() == before
    sheet.write_text("sample,stage\ninvented-a,E9\n")
    assert main([str(sheet), "--task", "t1"]) == 0


def test_cli_errors_hide_input_path(tmp_path, capsys):
    sheet = tmp_path / "missing.csv"
    assert main([str(sheet), "--task", "t1"]) == 2
    assert str(tmp_path) not in capsys.readouterr().err


def test_cli_explicit_atlas_header_mapping_and_delimiter(tmp_path, capsys):
    sheet = tmp_path / "invented-export.txt"
    sheet.write_text("specimen_key,assigned_day,Characteristics[age unit]\ninvented-a,8,days\n")
    assert main([str(sheet), "--task", "t2-embryo", "--delimiter", "comma",
                 "--sample-column", "specimen_key", "--stage-column", "assigned_day"]) == 0
    assert "OUTSIDE_BY_STAGE" in capsys.readouterr().out


def test_different_header_roles_cannot_select_the_same_column(tmp_path, capsys):
    sheet = tmp_path / "invented.csv"
    sheet.write_text("sample,stage\ninvented-a,E9\n")
    assert main([str(sheet), "--task", "t1", "--sample-column", "stage",
                 "--stage-column", "stage"]) == 2
    assert "distinct" in capsys.readouterr().err


def test_missing_sample_is_not_silently_merged():
    with pytest.raises(ValueError, match="identifier"):
        check_rows([{"sample": "NA", "stage": "E9"}], task="t1",
                   sample_column="sample", stage_column="stage")


def test_unknown_task_and_nonfinite_api_margin_rejected():
    with pytest.raises(ValueError, match="Unknown task"):
        check("E9", "typo")
    with pytest.raises(ValueError, match="finite"):
        check("E9", near=D("NaN"))
