# VE Community: stage checks and submission selection

`ve-stagecheck` checks a local public atlas sample sheet against the published
external-data stage windows. It groups technical records by sample, handles
explicit age intervals and uncertainty, and routes morphological or litter-based
ages to manual review. This repository contains code, documentation and invented
metadata fixtures only. It has no network client, challenge files, keys or
scientific outcomes.

Contents:

- `ve-stagecheck`: specimen-level stage-window triage for public atlas sample sheets (below).
- [tutorials/selection](tutorials/selection/README.md): synthetic-only tutorial on choosing what to submit after a large candidate search.
- [tutorials/proxy](tutorials/proxy/README.md): synthetic-only tutorial on whether a local proxy board is as hard as the test: horizon leakage, disjoint-step transfer and choosing between disagreeing proxies.

MIT licensed. Built for participants of the [Virtual Embryo Challenge](https://virtualembryo.ai/challenge). It uses no challenge data, and the challenge data stay under the challenge Terms of Use.

The [selection tutorial](tutorials/selection/README.md) uses only synthetic
scores to explain the winner's curse, sealed selection, paired comparisons and
expected improvement for a final best-of-k attempt. The [proxy tutorial](tutorials/proxy/README.md)
uses only simulated trajectories to show when a local proxy is easier than
the test. Both Python scripts run without runtime dependencies; their checks
run with the repository's pytest suite.

## Install and use

Python 3.10 or later; the tool has no third-party runtime dependencies.

```sh
python -m pip install -e '.[test]'
ve-stagecheck examples/synthetic_sample_sheet.csv --task t1
ve-stagecheck examples/synthetic.sdrf.txt --task t2-embryo
python -m pytest -q
```

For an actual public atlas sheet, supply exact headers when auto-detection is
ambiguous. `--sample-column` should identify the biological specimen, rather than
the sequencing run, so technical records are considered together.

```sh
ve-stagecheck public_sample_sheet.tsv --task t2-heart \
  --sample-column 'sample_id' --stage-column 'embryonic_day'
ve-stagecheck public.sdrf.txt --task t2-embryo \
  --sample-column 'Source Name' \
  --stage-column 'Characteristics[developmental stage]'
```

These are column-mapping examples, not promises that every atlas release uses
those headers. The Qiu `sample_sheet` and Pijuan-Sala SDRF must be inspected for
the appropriate specimen and age fields before use. Nothing is downloaded.
TSV and CSV work directly; other delimiters and Excel workbooks require a local
CSV/TSV export. UTF-8 BOMs, quoted fields and repeated SDRF source records work.

Output is TSV on stdout, with one decision per specimen, its source record
indices (starting at one after the header), reasons and rule references. Input
files are read only. Keep runtime output outside this repository; output sample
identifiers originate in the supplied metadata and are not anonymised.

| Decision | Meaning |
| --- | --- |
| `EXCLUDE` | A record clearly lies inside a protected stage window; exclude the whole specimen. |
| `REVIEW` | Resolve uncertain staging, conflicting records or condition questions before use. |
| `OUTSIDE_BY_STAGE` | Outside the encoded stage windows; other eligibility checks remain. |

Exit status: `0` means every specimen is outside by stage, `1` means at least one
exclusion or review, and `2` means a malformed sheet or invocation. The synthetic
examples deliberately contain exclusions, so their expected exit status is `1`.

## Staging and limitations

Accepted ages are explicit embryonic days, with optional `E`, `days` or `dpc`,
closed ranges such as `E9.4-E9.6`, and uncertainty such as `E9.5 +/- 0.1`.
Reversed ranges, approximate prose, unknown units, postnatal ages and missing
ages require review. Ages are never inferred from sample identifiers.

`--uncertainty-days` expands the stated interval on both sides. Use a value
justified by the source, not a guessed biological conversion. `--near-days`
adds an optional review band around boundaries; its default is zero, and this
user-selected band is **not** a published rule. Decimal arithmetic preserves
the published endpoint distinctions. An interval can touch allowed endpoints
while containing excluded ages in its interior.

Somite and Theiler columns with nonmissing values, including a somite count of
zero, require manual conversion even when a nominal E-day is also supplied.
Morphological labels in other staging columns and separately annotated non-day
age units also trigger review.
Litter-assigned ages require review; a litter identifier alone does not establish
how age was assigned. Preserve staging annotations in the sheet. Omitting an
uncertainty or staging-method column can hide information the checker needs.

Conflicting ages for one sample cause review, and any clear exclusion applies
to the entire sample. Disagreement between separate age columns also causes
review even when one column was explicitly selected. Task 3 always requires condition review: stage alone
cannot resolve allele equivalence, comparable stages or phenocopies. Non-WT
conditions and non-mouse species, when annotated, also trigger review. This
does not evaluate corpus contents, pretrained models, licensing, disclosure,
or organismal staging biology, and cannot certify external-data eligibility.
Ask the organisers about unresolved boundary cases before using a source.

## Sources and overlap

The [public rule notes](docs/PUBLIC_RULES.md) cite Rules §10 and the FAQ.
The [prior-art review](docs/PRIOR_ART.md) explains why the proposed count-recovery
and preflight tools were dropped and how this sheet-specific extension differs
from existing stage checkers. The [filing index](docs/FILINGS_INDEX.md) lists
existing public resources with filing dates and review status.

Credit: i-habib's [Data Safety Kit](https://github.com/i-habib/virtual-embryo-community/tree/main/data-safety)
addresses the underlying stage-rule problem; this is an independent implementation
of specimen-level sheet triage. The official [veckit](https://github.com/aristoteleo/veckit)
provides format/scoring infrastructure and is credited in the reuse notes;
this tool neither imports its scorer nor evaluates predictions.

Atlas references: Qiu et al. (2024), *A single-cell time-lapse of mouse prenatal
development from gastrula to birth*, [Nature](https://doi.org/10.1038/s41586-024-07069-w);
Pijuan-Sala et al. (2019), *A single-cell molecular map of mouse gastrulation and
early organogenesis*, [Nature](https://doi.org/10.1038/s41586-019-0933-9),
ArrayExpress accession E-MTAB-6967. No atlas records are bundled or used in tests.

Recheck the current rules and the live public register before using any source.
