# Is your proxy as hard as the test?

MIT licensed. For Virtual Embryo Challenge teams who score candidates on a
local proxy board built from released stages. Every trajectory, stage,
depth and score below is simulated. No challenge data or results from real
runs are used. Stages are integer times `0, 1, 2, ...`, not real ages.

A proxy is useful when it is as hard as the test in the ways that matter
for the candidates you compare. This lesson shows three ways a proxy can be
easier or different, a free check you can run before trusting an
extrapolation rule, and a way to choose between proxies that disagree.

Choosing among many candidates on one board is a separate problem. It is
covered by the [selection tutorial](../selection/README.md) in this
repository: winner's curse, sealed selection, paired comparisons.

## Run it

From the repository root, with Python 3.10 or later:

```sh
python tutorials/proxy/proxy.py
python -m pytest -q tests/test_proxy.py
```

The script uses only the standard library, reads no files, makes no network
requests and prints plain text. It runs in under a second. Fixed seeds make
the output repeatable. `--repeats` and `--genes` change Monte Carlo
precision, not the rules. With matplotlib installed,
`--plot /tmp/proxy.png` also writes a two-panel figure. Write it outside the
repository.

## 1. Horizon leakage

Each synthetic gene follows a smooth quadratic in time. Each stage is one
noisy specimen. Stages `0` to `4` are released. The test asks for stage `6`
from stage `4`: a horizon of two steps, with no observed stage in between.

Three boards:

| Board | Anchor to target | Stages a method can read | Observed after anchor |
| --- | --- | --- | --- |
| leaky proxy | 2 to 4 | 0, 1, 2, 3 | 3 |
| clean proxy | 2 to 4 | 0, 1, 2 | none |
| test | 4 to 6 | 0 to 4 | none |

The leaky proxy has the right nominal horizon. But a common way to build it
is to hold out stage `4` and keep every other released file. Stage `3`
stays readable, and it lies inside the target interval.

Two methods:

- `trend`: extrapolate linearly from the anchor and the stage before it.
  It never reads anything after the anchor.
- `reuse`: average the trend with the observed stage nearest to the target.
  It does not ask whether that stage would exist at test time.

Default output (mean squared error, lower is better):

```text
board        anchor->target  after anchor    trend    reuse  reuse-trend    sem
leaky proxy            2->4             3    1.340    0.757       -0.583  0.014
clean proxy            2->4          none    1.340    1.740        0.400  0.024
test                   4->6          none    1.342    1.929        0.587  0.027
```

On the leaky proxy, `reuse` looks much better, because its nearest stage is
one step from the target. On the test the nearest stage is the anchor, two
steps away, and `reuse` loses. The clean proxy gives the same sign as the
test. `trend` has the same error on both proxies because it never reads
stage `3`.

The leak is not about the held-out test stage. It is about information
inside the horizon that the test will not have. Any method that interpolates,
smooths over neighbours, or picks the closest reference can use it.

To build a clean proxy:

1. Set the proxy horizon equal to the test horizon, in steps or in time.
2. Take the latest released stage as the target and place the anchor one
   horizon earlier.
3. Remove every stage after the anchor from what the method can read.
   This includes files used for normalisation, gene selection and
   reference choices.

`clean_proxy(released, horizon)` does steps 2 and 3. `steps_after_anchor(proxy)`
returns the observed stages later than the anchor; an empty result is
required. [VEC SplitForge](https://github.com/random-guy-05/vec-splitforge)
(filed 2026-10-06) already builds past-to-future folds from a stage manifest
and audits them against the protected stage windows. Its extrapolation
folds use the immediately preceding stage as the anchor. This section adds
the horizon-length match and shows the cost of an observed stage inside the
horizon. Specimen leakage across a cell-level split is a different problem,
covered by [Reliability Lab 15](https://github.com/random-guy-05/vec-reliability-lab/blob/main/resources/en/15-specimen-holdout.md).

Background: rolling-origin evaluation for time series keeps every training
observation before the forecast origin and scores at a fixed horizon. See
Hyndman and Athanasopoulos,
[Forecasting: Principles and Practice, time series cross-validation](https://otexts.com/fpp3/tscv.html),
and Bergmeir and Benitez (2012),
[On the use of cross-validation for time series predictor evaluation](https://doi.org/10.1016/j.ins.2011.12.028),
Information Sciences 191:192-213.

## 2. Depth mismatch: covered elsewhere

Two stages sequenced at different depths differ in detection, so part of an
apparent stage-to-stage change can be depth. A proxy built from
depth-matched stages cannot show that effect. This topic is already covered
by [vec-t1-depth-covariate](https://github.com/drkeerthiharikrishnan-h/vec-t1-depth-covariate)
(filed 2026-10-29). It uses binomial thinning to a common depth to separate
the depth part of a step, builds a depth-structured local replica, and gives
a [depth-neutral fitting recipe](https://github.com/drkeerthiharikrishnan-h/vec-t1-depth-covariate/blob/main/depth_neutral_fit.md).
Thinning needs integer counts. For recovering counts from normalised
matrices, use i-habib's
[Count Recovery](https://github.com/i-habib/virtual-embryo-community/tree/main/count-recovery)
and check its assumptions. This tutorial does not repeat either.

Depth comes back in section 4 as one of the properties a proxy must match.

## 3. Disjoint-step transfer: a free falsification test

Many rules extrapolate the last observed step: `next = last + beta * (last - previous)`.
Such a rule assumes that the next step resembles the last one. You can test
that assumption on released stages alone, with no submission.

Fit on one released step and predict a later step that shares no stage
with it. With released stages `0` to `4`, compare step `1 -> 2` with step
`3 -> 4`, gene by gene. Their correlation estimates how much of a step
carries over to a later step.

The simulation gives each gene a persistent drift plus a new random
program at every step. `persistence` is the share of step variance carried
by the drift. It equals the true correlation between any two steps. Stage
`5` plays the test.

```text
persistence disjoint r adjacent r   beta    copy  extrap  shrunk
       0.80      0.726      0.686  0.703   1.061   0.646   0.551
       0.05      0.051      0.002  0.050   1.082   2.180   1.084
```

`copy` predicts stage `5` as stage `4`. `extrap` adds the full last step.
`shrunk` adds `beta` times the last step, where `beta` is the least-squares
slope of the later released step on the earlier one, clipped at zero.

With high persistence the disjoint correlation is high and extrapolation
helps. With low persistence the correlation is near zero, and adding the
full last step doubles the error relative to copying. The shrunk rule falls
back to copying. A weak disjoint correlation is a warning before any
submission is spent.

Use disjoint steps, not adjacent ones. Steps `2 -> 3` and `3 -> 4` share
stage `3`. Its specimen noise enters one step with a plus sign and the other
with a minus sign, so adjacent steps correlate lower. With specimen noise
SD `s`, unit step variance and persistence `p`, disjoint steps correlate at
`p / (1 + 2 s^2)` and adjacent steps at `(p - s^2) / (1 + 2 s^2)`. For pure
noise the adjacent value approaches `-0.5`. The test suite checks this case.

Limits: the check needs at least four released stages. It measures
gene-wise persistence, not cell-level structure. A high correlation on
released steps does not prove the unreleased step behaves the same. A low
one is strong evidence that the rule's premise fails on the data you have.
Steps of different lengths are not directly comparable.

## 4. When two proxies disagree

Candidate: extrapolate the last step (`beta = 1`). Baseline: copy the anchor.
The gain is baseline error minus candidate error, so positive means the
candidate helps. Three boards differ in three design properties, all known
before any scoring:

- horizon, in steps;
- the log depth change inside the last observed step. A stand-in shifts the
  shallower stage by a gene-specific amount. It is not a count model;
- specimens averaged per stage.

```text
board    horizon  depth specimens    gain     sd sign as test
proxy A        1   0.00         3   0.171  0.090         0.03
proxy B        2  -0.70         1  -0.449  0.287         0.95
test           2  -0.80         1  -0.611  0.394         0.97
Chosen before scoring, by design match only: proxy B
One factor changed from proxy A toward the test:
  A, test horizon    gain   0.838
  A, test depth      gain  -0.038
  A, test specimens  gain  -0.035
```

Proxy A says the candidate helps. Proxy B says it hurts. The test agrees
with B. `sign as test` is the share of repeats in which a board gives the
same sign as the mean test gain.

Decide which proxy to trust before looking at its scores.
`closer_design(test, proxies)` takes design descriptors only. It counts
matches on horizon, depth change and specimen count, and returns the proxy
with the most matches. On a tie or no match it returns `None`: neither proxy
should decide, and a matched proxy should be built. The function signature
has no score argument, and the tests check that.

The ablation shows why all three properties matter together. A longer
horizon alone makes extrapolation look better in this model, because the
persistent drift adds up. At horizon one, a depth step or a single specimen
alone only removes the small gain. The candidate adds `h` times the last
step, so a depth shift in that step enters `h` times and anchor noise
enters `1 + h` times. Their costs grow with `h^2`. Only the board that
matches all three properties shows the net loss.

Do not pick the proxy whose answer you prefer, the one with the larger
gain, or the one with the smaller spread. Proxy A has the smallest spread
because it averages three specimens. It is precise about a different
problem.

Related work. [Reliability Lab 14](https://github.com/random-guy-05/vec-reliability-lab/blob/main/resources/en/14-proxy-rank-reversal.md)
(filed 2026-10-21) shows that a local board can reverse official rankings
and how to measure that with tie-aware rank correlation.
[Val-to-test](https://github.com/i-habib/virtual-embryo-community/tree/main/val-to-test)
(filed 2026-10-15) measures how a validation-like ranking carries over to a
test-like pseudo pair with a longer step.
[VecBench](https://github.com/AJK0921/vecbench) (filed 2026-08-29) measures
whether a local holdout can rank candidates at all, and
[vec-board-calibration](https://github.com/drkeerthiharikrishnan-h/vec-board-calibration)
(filed 2026-10-29) compares a local replica with returned board values.
This section adds a design-only rule for choosing between proxies before
scoring.

## 5. Checklist

1. **Match the horizon.** Same number of steps, or the same time gap, from
   anchor to target as the test.
2. **Keep observed steps out of the target interval.** Nothing after the
   anchor may be read by the method, including preprocessing.
3. **Match or model depth.** If the test's last observed step crosses a
   depth change, the proxy's should too, or remove depth first. See
   [vec-t1-depth-covariate](https://github.com/drkeerthiharikrishnan-h/vec-t1-depth-covariate).
4. **Count specimens.** One specimen per stage limits what any proxy can
   say. Specimen noise enters every step, and extrapolation amplifies it.
   Averaging specimens on a proxy when the test has one makes the proxy
   easier, not better. Proxy differences smaller than their repeat spread
   are not evidence.
5. **Run the disjoint-step check** before trusting a rule that extrapolates
   the last step.
6. **Freeze the proxy before scoring candidates.** Choose the design, the
   stages and the metric first. Changing the proxy after seeing results
   is selection on the proxy itself; see Cawley and Talbot (2010),
   [On Over-fitting in Model Selection and Subsequent Selection Bias in Performance Evaluation](https://jmlr.org/papers/v11/cawley10a.html),
   JMLR 11:2079-2107, and the [selection tutorial](../selection/README.md).

## Limits

The trajectories are smooth, gene-wise and independent across genes. Real
expression has cell-type composition, correlated genes and nonlinear
programs. The depth term is a fixed gene-wise shift, not a sampling model.
Mean squared error stands in for the challenge metrics, which can respond
differently. The examples show mechanisms. They do not predict any board
score or the size of any effect on real data.
