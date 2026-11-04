"""Conservative sample-sheet triage. Public rules snapshot: 2026-11-04.

Independent implementation extending the problem addressed by i-habib's
Data Safety Kit. See docs/PRIOR_ART.md and docs/PUBLIC_RULES.md.
"""
from __future__ import annotations

import argparse
import csv
import io
import re
import sys
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from pathlib import Path

RULES_DATE = "2026-11-04"
RULES_URL = "https://virtualembryo.ai/challenge/rules"
FAQ_URL = "https://virtualembryo.ai/challenge/faq"
RULE_REFS = (
    f"{RULES_URL} s10 L89-L102; {FAQ_URL} L51-L57 "
    "(rendered text, 2026-11-04)"
)
TASKS = ("t1", "t2-heart", "t2-embryo", "t3")
D = Decimal
ZERO = D(0)


@dataclass(frozen=True)
class Window:
    lo: Decimal
    hi: Decimal
    closed_lo: bool = False
    closed_hi: bool = False

    def contains(self, value: Decimal) -> bool:
        return (value > self.lo or (self.closed_lo and value == self.lo)) and (
            value < self.hi or (self.closed_hi and value == self.hi)
        )

    def intersects(self, lo: Decimal, hi: Decimal) -> bool:
        left, right = max(lo, self.lo), min(hi, self.hi)
        return left < right or (left == right and self.contains(left))


# Public rules constants, not measurements. FAQ explicitly covers heart extrapolation.
EXTRAPOLATION = Window(D("9.5"), D("13.5"), closed_hi=True)
WINDOWS = {
    "t1": (EXTRAPOLATION,),
    "t2-heart": (Window(D("8.25"), D("8.75")), EXTRAPOLATION),
    "t2-embryo": (Window(D("7.25"), D("8.0")),),
    "t3": (),
}


@dataclass(frozen=True)
class Age:
    lo: Decimal
    hi: Decimal


@dataclass(frozen=True)
class Verdict:
    sample: str
    rows: tuple[int, ...]
    decision: str
    nominal_relation: str
    reasons: tuple[str, ...]


NUMBER = r"(?:\d+(?:\.\d+)?|\.\d+)"
TOKEN = rf"[eE]?\s*({NUMBER})"
UNIT = r"(?:\s*(?:dpc|days?(?:\s+post\s+coitum)?|embryonic\s+days?))?"
SINGLE = re.compile(rf"{TOKEN}{UNIT}", re.IGNORECASE)
RANGE = re.compile(rf"{TOKEN}\s*(?:-|–|—|to)\s*{TOKEN}{UNIT}", re.IGNORECASE)
PLUSMINUS = re.compile(rf"{TOKEN}\s*(?:\+/-|±)\s*({NUMBER}){UNIT}", re.IGNORECASE)
MISSING = {"", "na", "n/a", "nan", "none", "null", "unknown", "not available", "-"}


def present(value: str) -> bool:
    return value.strip().lower() not in MISSING


def parse_age(text: str) -> Age | None:
    """Full-string E-day grammar only; never infer an age from an identifier."""
    text = text.strip()
    match = PLUSMINUS.fullmatch(text)
    if match:
        center, delta = map(D, match.groups())
        lo, hi = center - delta, center + delta
    else:
        match = RANGE.fullmatch(text)
        if match:
            lo, hi = map(D, match.groups())
        else:
            match = SINGLE.fullmatch(text)
            if not match:
                return None
            lo = hi = D(match.group(1))
    # Mouse embryonic days only. Postnatal stages need an explicit external mapping.
    if not (D(0) <= lo <= hi <= D(30)):
        return None
    return Age(lo, hi)


def relation(age: Age, task: str) -> str:
    windows = WINDOWS[task]
    if any(w.contains(age.lo) and w.contains(age.hi) for w in windows):
        return "INSIDE"
    if any(w.intersects(age.lo, age.hi) for w in windows):
        return "CROSSES"
    return "OUTSIDE"


def normalize_header(value: str) -> str:
    return re.sub(r"[\s_]+", " ", value.strip().lower())


ALIASES = {
    "sample": {"sample", "sample id", "sample name", "source name", "embryo", "embryo id"},
    "stage": {"stage", "age", "embryonic day", "developmental stage", "day",
              "characteristics[age]", "characteristics[developmental stage]"},
}


def choose_column(headers: Sequence[str], role: str, explicit: str | None) -> str:
    if explicit is not None:
        if explicit not in headers:
            raise ValueError(f"Selected {role} column is absent; inspect the sheet header.")
        return explicit
    candidates = [h for h in headers if normalize_header(h) in ALIASES[role]]
    if len(candidates) != 1:
        raise ValueError(f"Select --{role}-column explicitly: header inference is missing or ambiguous.")
    return candidates[0]


def read_sheet(path: Path, delimiter: str | None = None) -> tuple[list[str], list[dict[str, str]]]:
    """Read local CSV/TSV/SDRF metadata, preserving repeated specimen rows."""
    with path.open(encoding="utf-8-sig", newline="") as handle:
        first = handle.readline()
        handle.seek(0)
        if delimiter is None:
            if path.suffix.lower() == ".csv":
                delimiter = ","
            elif "\t" in first:
                delimiter = "\t"
            else:
                delimiter = ","
        reader = csv.DictReader(handle, delimiter=delimiter, strict=True)
        headers = reader.fieldnames or []
        if not headers or any(not h.strip() for h in headers):
            raise ValueError("The sheet needs a nonempty header for every column.")
        normalized = [normalize_header(h) for h in headers]
        if len(set(normalized)) != len(normalized):
            raise ValueError("Duplicate headers are ambiguous; select or rename columns first.")
        rows = []
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError("Ragged sheet: every record must match the header width.")
            rows.append(row)
        if not rows:
            raise ValueError("The sheet has no sample rows.")
        return headers, rows


def annotation_flags(row: Mapping[str, str], stage_column: str) -> list[str]:
    flags = []
    selected_age = parse_age(row[stage_column])
    for header, value in row.items():
        if not present(value):
            continue
        h, v = normalize_header(header), value.strip().lower()
        # A litter identifier alone is not a claim that the age was assigned by litter.
        staging_field = header == stage_column or any(term in h for term in ("stage", "age", "staging"))
        if header != stage_column and h in ALIASES["stage"]:
            alternate_age = parse_age(value)
            if alternate_age is not None and selected_age is not None and alternate_age != selected_age:
                flags.append("Alternative age columns disagree; reconcile the source annotations.")
        if "somite" in h or "theiler" in h or (
            staging_field and re.search(r"somite|theiler|\bts\s*\d", v)
        ):
            flags.append("Alternative morphological staging needs manual conversion.")
        if "litter age" in h or "litter stage" in h or (
            (staging_field or any(word in h for word in ("basis", "method", "age source")))
            and "litter" in v
        ):
            flags.append("Litter-assigned age needs specimen-level review.")
        if staging_field and "unit" in h and v not in {
            "day", "days", "d", "dpc", "embryonic day", "embryonic days", "days post coitum"
        }:
            flags.append("Separate age-unit annotation is not embryonic days; convert manually.")
        if (any(word in h for word in ("genotype", "perturbation", "phenocopy", "allele"))
                and v not in {"wt", "wild type", "wild-type", "control", "mus musculus"}):
            flags.append("Perturbation or allele requires a separate condition review.")
        if (h in {"species", "organism", "characteristics[organism]"}
                and v not in {"mus musculus", "mouse"}):
            flags.append("Non-mouse or unresolved species needs manual review.")
    return list(dict.fromkeys(flags))


def inspect_row(row: Mapping[str, str], task: str, stage_column: str,
                uncertainty: Decimal, near: Decimal) -> tuple[str, str, list[str], Age | None]:
    flags = annotation_flags(row, stage_column)
    age = parse_age(row[stage_column])
    if age is None:
        return "REVIEW", "UNRESOLVED", flags + ["No unambiguous embryonic-day age."], None
    age = Age(max(D(0), age.lo - uncertainty), age.hi + uncertainty)
    if task == "t3":
        # Stage alone cannot certify allele equivalence, comparability or phenocopy.
        return "REVIEW", "CONDITION_DEPENDENT", flags + [
            "Task 3 needs genotype, comparable-stage and phenocopy review."
        ], age
    nominal = relation(age, task)
    if near and any(min(abs(age.lo - edge), abs(age.hi - edge)) <= near
                    for w in WINDOWS[task] for edge in (w.lo, w.hi)):
        flags.append("Within the user-selected boundary review band; this band is not a rule.")
    if nominal == "CROSSES":
        flags.append("The age interval crosses a protected boundary; resolve before use.")
    if flags:
        return "REVIEW", nominal, flags, age
    if nominal == "INSIDE":
        return "EXCLUDE", nominal, ["Inside a published protected stage window."], age
    return "OUTSIDE_BY_STAGE", nominal, [
        "Outside the stage windows only; licence, provenance and disclosure still apply."
    ], age


def check_rows(rows: Iterable[Mapping[str, str]], *, task: str, sample_column: str,
               stage_column: str, uncertainty: Decimal = ZERO, near: Decimal = ZERO
               ) -> list[Verdict]:
    if task not in TASKS:
        raise ValueError("Unknown task; select a supported task name.")
    if not all(v.is_finite() and v >= 0 for v in (uncertainty, near)):
        raise ValueError("Review margins must be finite and nonnegative.")
    grouped: dict[str, list[tuple]] = {}
    for index, row in enumerate(rows, start=1):
        sample = row[sample_column].strip()
        if not present(sample):
            raise ValueError("Every row needs a sample identifier; missing identifiers cannot be grouped.")
        item = inspect_row(row, task, stage_column, uncertainty, near)
        grouped.setdefault(sample, []).append((index, *item))
    if not grouped:
        raise ValueError("The sheet has no sample rows.")
    verdicts = []
    for sample, items in grouped.items():
        decisions = {item[1] for item in items}
        ages = {item[4] for item in items}
        reasons = list(dict.fromkeys(reason for item in items for reason in item[3]))
        if len(ages) > 1:
            reasons.append("Conflicting age records for one sample; review the entire sample.")
        # A clearly prohibited record excludes the entire specimen, even if other records differ.
        if "EXCLUDE" in decisions:
            decision = "EXCLUDE"
        elif "REVIEW" in decisions or len(ages) > 1:
            decision = "REVIEW"
        else:
            decision = "OUTSIDE_BY_STAGE"
        relations = {item[2] for item in items}
        verdicts.append(Verdict(sample, tuple(item[0] for item in items), decision,
                                next(iter(relations)) if len(relations) == 1 else "MIXED",
                                tuple(reasons)))
    return verdicts


def render(verdicts: Sequence[Verdict], task: str) -> str:
    buf = io.StringIO(newline="")
    writer = csv.writer(buf, delimiter="\t", lineterminator="\n")
    writer.writerow(("sample", "records", "task", "decision", "nominal_relation", "reason", "source"))
    for v in verdicts:
        # Prevent spreadsheet formula evaluation when a public sample ID starts with an operator.
        sample = "'" + v.sample if v.sample.startswith(("=", "+", "-", "@")) else v.sample
        writer.writerow((sample, ",".join(map(str, v.rows)), task, v.decision,
                         v.nominal_relation, " ".join(v.reasons), RULE_REFS))
    return buf.getvalue()


def margin(value: str) -> Decimal:
    try:
        number = D(value)
    except InvalidOperation as error:
        raise argparse.ArgumentTypeError("Margin must be a finite nonnegative day value.") from error
    if not number.is_finite() or number < 0:
        raise argparse.ArgumentTypeError("Margin must be a finite nonnegative day value.")
    return number


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sheet", type=Path, help="Local public atlas CSV, TSV or SDRF")
    parser.add_argument("--task", required=True, choices=TASKS)
    parser.add_argument("--sample-column", help="Exact sample/specimen header; auto if unambiguous")
    parser.add_argument("--stage-column", help="Exact embryonic-day header; auto if unambiguous")
    parser.add_argument("--delimiter", choices=("comma", "tab"))
    parser.add_argument("--uncertainty-days", type=margin, default=D(0),
                        help="Expand every nominal age by this documented uncertainty")
    parser.add_argument("--near-days", type=margin, default=D(0),
                        help="Optional conservative review band; not a challenge rule")
    args = parser.parse_args(argv)
    try:
        delimiter = {"comma": ",", "tab": "\t"}.get(args.delimiter)
        headers, rows = read_sheet(args.sheet, delimiter)
        sample = choose_column(headers, "sample", args.sample_column)
        stage = choose_column(headers, "stage", args.stage_column)
        if sample == stage:
            raise ValueError("Sample and stage must use distinct columns.")
        verdicts = check_rows(rows, task=args.task, sample_column=sample, stage_column=stage,
                              uncertainty=args.uncertainty_days, near=args.near_days)
    except (OSError, UnicodeError):
        print("Cannot read the local metadata sheet; check access and UTF-8 encoding.", file=sys.stderr)
        return 2
    except (ValueError, csv.Error) as error:
        print(str(error), file=sys.stderr)
        return 2
    sys.stdout.write(render(verdicts, args.task))
    return int(any(v.decision != "OUTSIDE_BY_STAGE" for v in verdicts))


if __name__ == "__main__":
    raise SystemExit(main())
