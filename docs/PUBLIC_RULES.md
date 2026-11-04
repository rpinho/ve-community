# Public rule snapshot for ve-stagecheck

Reviewed 2026-11-04. Stage numbers here are public rule constants, not observed
values. Line references identify the web tool's rendered-text view on that date;
they are citation aids, not permanent website line numbers. The tool is offline
and does not silently update this snapshot.

| Task | Protected external-data stage windows | Source lines |
| --- | --- | --- |
| `t1` | `(E9.5, E13.5]` | Rules §10 L94–L96; FAQ L54–L55 |
| `t2-heart` | `(E8.25, E8.75)` and `(E9.5, E13.5]` | Rules §10 L94–L98; FAQ L54–L56 explicitly includes heart extrapolation |
| `t2-embryo` | `(E7.25, E8.0)` | Rules §10 L97–L99; FAQ L56 |
| `t3` | Condition-dependent; no stage-only clearance | Rules §10 L91, L100; FAQ L52, L57 |

Sources: [Rules §10](https://virtualembryo.ai/challenge/rules),
[FAQ, external public data](https://virtualembryo.ai/challenge/faq).
Explicit extrapolation and interpolation windows govern this implementation;
do not recalculate them from the earlier general midpoint wording. Round brackets
exclude the endpoint; square brackets include it in the prohibited interval.

Morphological staging must be assessed on the same developmental scale; an E-day
label alone is insufficient (Rules L92, L98, L102; FAQ L53). Mixed external
resources require removal of protected samples before training (Rules L101).
Disclosure remains required (Rules L89); pretrained exposure and comparable
alleles or phenocopies require separate scrutiny (Rules L90, L100).

Official validation access is a specific P3 exception, **not** permission to use
external samples from these windows (Rules L86–L87; FAQ L51). This tool handles
public external metadata only and implements no validation-data exception.
