# Overlap review and scope decision

Reviewed on 2026-11-04. The [award page](https://virtualembryo.ai/challenge/community)
asks contributors to check earlier filings. The [resources page](https://virtualembryo.ai/challenge/resources)
loads its entries client-side from a public read-only
[register](https://api.virtualembryo.ai/challenge/community/register).
The static HTML alone appears empty. Both reviewed and pending entries were
examined; [the link index](FILINGS_INDEX.md) retains their public filing metadata.
Listing and acceptance do not establish correctness or broad eligibility.

## Existing work relevant to the proposals

| Existing work | Public source | Scope and decision |
| --- | --- | --- |
| Count Recovery | [i-habib/count-recovery](https://github.com/i-habib/virtual-embryo-community/tree/main/count-recovery), filed 2026-10-09 | Integer recovery, library depths and round-trip checks already exist. Drop `ve-counts`. |
| Depth covariate method | [vec-t1-depth-covariate](https://github.com/drkeerthiharikrishnan-h/vec-t1-depth-covariate), filed 2026-10-29 | Prior source for the depth-confounding idea and depth-neutral fitting; no outcome claims or artifacts reused. |
| Reliability Lab | [count-recovery lesson](https://github.com/random-guy-05/vec-reliability-lab/blob/main/resources/en/01-count-recovery.md) | Synthetic count-identifiability lessons already exist; drop the proposed replacement on that topic too. |
| Submission Validator | [xxx12e/vec_submit_check](https://github.com/xxx12e/vec-community-kit/tree/main/vec_submit_check), filed 2026-09-06 | Board panels, structure and local format checks already exist. Drop `ve-preflight`. |
| Evidence Skeleton | [xxx12e/vec_agent_evidence](https://github.com/xxx12e/vec-community-kit/tree/main/vec_agent_evidence), filed 2026-09-06 | Evidence collection, configuration locks and packaging already exist. |
| Evidence Check | [vec-evidence-check](https://github.com/bestdeeplearning-star/vec-evidence-check), filed 2026-09-14 | Offline evidence audits already exist, including a later per-submission quota correction. |
| Data Safety Kit | [i-habib/data-safety](https://github.com/i-habib/virtual-embryo-community/tree/main/data-safety), filed 2026-10-09 | Checks individual stages/ranges and writes disclosures. Generic stage checks would duplicate it. |
| External eligibility tool | [participation workbench](https://random-guy-05.github.io/vec-participation-workbench/tools/external-data-eligibility/), filed 2026-10-15 | Existing interactive eligibility guidance; no claim to invent the stage-checking idea. |
| Atlas catalogs | [External Data Catalog](https://github.com/i-habib/external-data-catalog) and [EmbryoMatch](https://github.com/ShaoliZhao/virtual-embryo-EmbryoMatch) | Existing source registries and preparation tools; do not create another catalog. |
| Normalisation utilities | [Panel Normalise](https://github.com/i-habib/virtual-embryo-community/tree/main/panel-normalise) | Existing raw-count/panel transformations. |

Public source interfaces were inspected at
`i-habib/virtual-embryo-community@d5b08ee31e42448b363b7add43267819b29aef2a`,
`xxx12e/vec-community-kit@0569db3233dd0729d8dd44adc18e5f33ca00436a`,
`vec-t1-depth-covariate@750955618e6930c2aa22852ef9bcc0bafc36d81c`, and
`vec-reliability-lab@62186e8214b54ce136f6553df769481f863fdb81`.
Code was read, not executed or copied. No scientific files or outcome tables
were fetched. The existing public-repository findings supplied for this job
were also used as a deduplication input; none of their private analysis is retained.

Other known repositories include tubager/vec-t2, wqty123/virtual-embryo-challenge-methods,
ShaoliZhao/VirtualEmbryoR, VictorRagazzi/virtual-embryo, and
Boom5426/Virtual_Embryo_Challenge_EigenLab (baselines or modelling methods);
random-guy-05's navigator, pathway builder, workshops and rules trainer
(education/navigation); aaronredzsir-star/perturbalign (schema adaptation);
and aristoteleo/veckit plus its public forks (scorer infrastructure). Those
families are outside the retained tool's specimen-level metadata workflow.
This review establishes bounded overlap, not proof of global novelty.

## Retained extension

`ve-stagecheck` adds a batch workflow for local CSV/TSV/SDRF sample sheets:
explicit header mapping, repeated technical-record grouping by biological sample,
interval uncertainty, a whole-sample exclusion rule, and conflicting-age review.
It preserves somite/Theiler and litter-assignment warnings rather than turning
nominal ages into automatic clearance. The current FAQ also makes the heart
extrapolation window explicit; the inspected older Data Safety Kit source
deferred those ages to organisers.

This is an extension of existing external-data safety work. It does not claim
that rule-window logic itself is new, and award eligibility remains for the
organisers to decide. One focused tool is retained rather than filing copies
of count recovery or preflight tooling. Consider contributing the adapters
upstream after maintainer review.

## General depth caveat

Library depth can change detection rates and apparent expression distributions,
so differences in sequencing depth can confound stage-to-stage change even after
normalisation. Neither apparent shifts nor depth equalisation establish a
biological mechanism. The depth idea is credited to the linked depth-covariate
repository; no observations from any challenge run are included here.

The minimum-positive-value recovery premise is conditional: a row must actually
contain a one-molecule observation. Integer ratios alone do not prove that
condition; multiplying a count vector by a common integer preserves its
normalised profile. For a filtered gene panel, recovered panel counts also need
not sum to the original full-library total. Use an established recovery tool
with its assumptions checked rather than claim unconditional reconstruction.
