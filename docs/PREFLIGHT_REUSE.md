# Reuse existing preflight tooling

The proposed `ve-preflight` CLI was dropped because earlier filings already cover
it. Start with the credited [Submission Validator](https://github.com/xxx12e/vec-community-kit/tree/main/vec_submit_check)
and [Evidence Skeleton](https://github.com/xxx12e/vec-community-kit/tree/main/vec_agent_evidence),
or [Evidence Check](https://github.com/bestdeeplearning-star/vec-evidence-check).
Check current public rules before use.

The official [veckit loader](https://github.com/aristoteleo/veckit/blob/cc3859ab2a2305efc4544815d3f44f660efa77d4/score_h5ad.py#L54)
checks panel order, finite nonnegative expression, and coordinate structure.
Its `_prediction` code (lines 54–89 at that commit) is useful source credit,
but a local scorer is not a complete portal-contract validator. In particular,
allowing a prediction to define its own panel cannot verify that panel against
the official contract. Use the existing validator and current published panel;
do not invoke a scorer with held-out reference files as a format check.

For Agent Team evidence, check these items against
[Rules §9](https://virtualembryo.ai/challenge/rules), rendered lines L72–L80,
and prize requirements at L137:

- The trajectory plus a second distinct kind is the scoring minimum. Multiple
  files of one kind do not supply a second kind.
- File all trajectory, prompts and harness material for prize verification.
- Each file is at most 200 MiB, and their sum at most 2 GiB **per submission**.
  Evidence for another submission does not reduce this allowance.
- Record the actual framework and model. Evidence must describe the run that
  produced the prediction and preserve its configuration lock.

The FAQ's shorter Agent Team answer has older minimum-evidence wording;
use the explicit current Rules §9 requirement above and resolve discrepancies
with the organisers.

For an eventual **opt-in remote format-only request**, omit the
`Idempotency-Key` header entirely, and explicitly set `format_only=true`.
Do not send an empty header or reuse a scored-submission idempotency key.
The [public OpenAPI schema](https://api.virtualembryo.ai/openapi.json) identifies
the format flag and optional header on `POST /challenge/submissions`; it does
not document their interaction. The omission is an operational precaution,
not a guarantee inferred from that schema. No remote request was tested here.
Rules §11 L115–L118 describe format checks and the separate prediction size cap.

This repository implements no upload or API path and reads no credentials.
An eventual uploader must require explicit opt-in, obtain a key from an
environment variable and verify current portal behaviour before use.
