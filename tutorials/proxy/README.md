# Is your proxy as hard as the test?

MIT licensed. For Virtual Embryo Challenge teams who score candidates on a
local proxy board built from released stages. This is guidance for
temporal extrapolation: predicting a later stage from earlier ones, with
nothing observed after the anchor. The challenge also has
[interpolation tasks](https://virtualembryo.ai/challenge/tasks), where later
bracketing observations are legitimate inputs; a proxy for those must keep
the permitted bracketing geometry instead. Every trajectory, stage,
depth and score below is simulated. No challenge data or results from real
runs are used. Stages are integer times `0, 1, 2, ...`, not real ages.

A proxy is useful when it is as hard as the test in the ways that matter
for the candidates you compare. This lesson shows three ways a proxy can be
easier or different, how to check a rule that extrapolates the last step,
and a conservative screen for choosing between proxies that disagree.

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
requests and prints plain text. It runs in about a second. Fixed seeds make
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

`clean_proxy(released, horizon)` returns the stage set for steps 2 and 3;
enforcing it in preprocessing is the team's job. `steps_after_anchor(proxy)`
returns the observed stages later than the anchor; an empty result is
required. [VEC SplitForge](https://github.com/random-guy-05/vec-splitforge)
(filed 2026-10-06) already builds past-to-future folds from a stage manifest
and audits them against the protected stage windows. Its extrapolation
folds use the immediately preceding stage as the anchor. Arnav Mana's
[VEC Chronological Holdout Split Auditor](https://random-guy-05.github.io/vec-participation-workbench/tools/temporal-split-audit/)
(filed 2026-10-15) checks that validation and test stages in a split
manifest are strictly later than training stages. This section adds
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

## 3. Checking a last-step rule

Many rules extrapolate the last observed step: `next = last + beta * (last - previous)`.
Such a rule assumes that the next step resembles the last one. The main
check is to run the actual rule, at the intended horizon, as a rolling-origin
forecast on released stages: for each origin, fit on stages up to it,
predict the stage one horizon later, and compare with copying the origin.
Two step correlations help explain the result.

- Disjoint steps. With released stages `0` to `4`, correlate step `1 -> 2`
  with step `3 -> 4`, gene by gene. The steps share no stage, so specimen
  noise does not enter the correlation. It is a diagnostic of longer-lag
  persistence. It speaks to the next step only if dependence is equal
  across lags.
- Adjacent steps. Correlate step `2 -> 3` with step `3 -> 4`. A last-step
  rule is itself adjacent: it predicts one step from the step right before
  it. So its coefficient must come from an adjacent pair. `beta` below is
  the least-squares slope of step `3 -> 4` on step `2 -> 3`, clipped at zero.

The simulation gives each gene a persistent drift plus a new random
program at every step. `persistence` `p` is the share of step variance
carried by the drift. It equals the true correlation between any two steps,
at every lag. Specimen noise has SD `s`. Stage `5` plays the test. The last
row is a counterexample in which dependence is not equal across lags: each
step is `A cos(k pi/4) + B sin(k pi/4)` with independent per-gene `A` and
`B`, so true steps one apart correlate at 0.71 and steps two apart at 0.

```text
case            disjoint r adjacent r   beta   copy extrap shrunk    gain     sd released
p=0.80 s=0.2         0.737      0.705  0.705  1.073  0.643  0.550   0.524  0.050    0.526
p=0.05 s=0.2         0.044      0.015  0.026  1.079  2.196  1.081  -0.002  0.004   -0.001
p=0.30 s=0.6         0.167     -0.034  0.001  1.707  3.493  1.707  -0.000  0.000   -0.001
rotating s=0.2      -0.001      0.622  0.625  1.080  0.816  0.661   0.418  0.056    0.411
```

`copy` predicts stage `5` as stage `4`. `extrap` adds the full last step.
`shrunk` adds `beta` times the last step. Errors are mean squared errors
on stage `5`. `gain` is `copy` minus `shrunk`, and `sd` is its standard
deviation across the 20 repeats of 800 genes. `released` is the
rolling-origin gain of the shrunk rule on released stages only (origins
`2` and `3`), with no test stage read.

With unit step variance, the adjacent slope estimates
`(p - s^2) / (1 + 2 s^2)` and the disjoint correlation estimates
`p / (1 + 2 s^2)`. The shared stage's noise enters adjacent steps with
opposite signs, and it enters the rule's own error the same way. So
extrapolation can help only when persistence exceeds specimen-noise
variance, `p > s^2`. The third row shows the gap: the disjoint correlation
is positive, but `p = 0.3` is below `s^2 = 0.36`. The adjacent slope clips
to zero, and the shrunk rule falls back to copying. Shrinking by the
disjoint value would lose to copying. For pure noise the adjacent
correlation approaches `-0.5`. The test suite checks these cases.

The rotating row shows the other gap. The disjoint correlation is near
zero, yet the shrunk rule beats copying. A disjoint check alone would have
rejected a rule that helps. In every row the rolling-origin gain on
released stages is close to the gain on the test stage.

Limits: these checks need at least four released stages. They measure
gene-wise behaviour, not cell-level structure. Good results on released
steps do not prove the unreleased step behaves the same. Steps of
different lengths are not directly comparable. With few released stages
the rolling-origin check has few origins, so report its spread.

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
Design screen, before scoring: proxy B
One factor changed from proxy A toward the test:
  A, test horizon    gain   0.838  sem 0.048
  A, test depth      gain  -0.038  sem 0.018
  A, test specimens  gain  -0.035  sem 0.016
Counterexample for counting matches (same test):
  proxy C  horizon 2 depth  0.00 specimens 1  matches 2  gain   0.142  sem 0.050
  proxy D  horizon 1 depth -0.80 specimens 2  matches 1  gain  -0.082  sem 0.019
  counting matches picks: proxy C
  design screen picks:    abstain
```

Proxy A says the candidate helps. Proxy B says it hurts. The test agrees
with B. `sign as test` is the share of repeats in which a board gives the
same sign as the mean test gain.

Decide which proxy to trust before looking at its scores.
`design_screen(test, proxies)` is a conservative design screen. It takes
design descriptors only. Declare three criteria before scoring: the same
horizon, a depth change within a tolerance of the test's (0.15 log units
here, an illustrative choice), and the same specimen count. The screen
returns a proxy only if exactly one proxy meets all three. Otherwise it
returns `None` and abstains: no proxy settles the question, and a matched
proxy should be built. The function signature has no score argument, and
the tests check that. Meeting all three criteria does not guarantee that a
proxy transfers to the test; it only rules out the mismatches named here.

Counting matches is not a safe substitute. In the counterexample, proxy C
meets two criteria and proxy D one, so a most-matches rule picks C. C's
mean gain is positive and the test's is negative. The screen abstains.

The ablation shows how the properties interact. A longer horizon alone
makes extrapolation look better in this model, because the persistent
drift adds up. The candidate adds `h` times the last step, so a depth shift
in that step enters `h` times and anchor noise enters `1 + h` times. Their
costs grow with `h^2`. The combined mismatch produces a larger loss in this
example; depth-only and specimen-only changes from proxy A also produce
small negative mean gains.

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
This section adds a design-only screen for choosing between proxies before
scoring, which abstains unless one proxy matches the test on every
predeclared criterion.

## 5. Checklist

1. **Match the horizon.** Same number of steps, or the same time gap, from
   anchor to target as the test.
2. **Keep observed steps out of the target interval.** For an extrapolation
   test with no post-anchor observations, nothing after the anchor may be
   read by the method, including preprocessing. For an interpolation test,
   keep the permitted bracketing geometry instead.
3. **Match or model depth.** If the test's last observed step crosses a
   depth change, the proxy's should too, or remove depth first. See
   [vec-t1-depth-covariate](https://github.com/drkeerthiharikrishnan-h/vec-t1-depth-covariate).
4. **Count specimens.** One specimen per stage limits what any proxy can
   say. Specimen noise enters every step, and extrapolation amplifies it.
   Averaging specimens on a proxy when the test has one makes the proxy
   easier, not better. Report uncertainty for paired differences across
   independent evaluation units. Distinguish uncertainty in mean
   performance from variability on a future specimen; repeated seeds do
   not create additional biological specimens. See
   [paired seeds](../selection/README.md#4-paired-seeds) in the selection
   tutorial.
5. **Check a last-step rule by rolling-origin forecasts** of the actual
   rule, at the intended horizon, on released stages. Take its coefficient
   from an adjacent pair of steps, clipped at zero. Read the disjoint-step
   correlation as a persistence diagnostic, not a verdict.
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
