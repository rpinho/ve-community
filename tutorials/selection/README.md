# Choosing what to submit without fooling yourself

MIT licensed. For Virtual Embryo Challenge teams, especially Agent
Team entries that generate many candidates. Every quality, score, panel and
random stream below is invented. No challenge data or results from real runs
are used. Scores have arbitrary higher-is-better units.

This lesson separates three questions. Which candidate looks best on an
exploratory board? Is there evidence of improvement on an independent proxy?
Which final attempt increases an already protected best score? Those questions
can have different answers.

## Run it

From the repository root, with Python 3.10 or later:

```sh
python tutorials/selection/selection.py
python -m pytest -q
```

Install the test extra first if pytest is absent: `python -m pip install -e '.[test]'`.
The script uses only the standard library. It prints text, reads no input files,
makes no network requests and writes no results. Default work is a few million
scalar Gaussian draws. Fixed seeds make the examples repeatable. For a shorter
smoke run, use `--repeats 30 --trials 100`. These options change Monte Carlo
precision, not the decision rules. There is no parameter search.

## 1. The winner's curse

Give every candidate true quality `q = 0.50`. A small board has 16 independent
units, with independent Gaussian noise of SD `0.20` per candidate and unit.
Its mean score has noise SD `s = 0.20 / sqrt(16) = 0.05`.
The code draws from this exact distribution of board means; it does not need
to materialize all 16 units. Each trial scores all N candidates on one board
and selects its maximum. Repeating the whole trial estimates selection bias.

For independent, equal-quality Gaussian candidates, the leading extreme-value
approximation is

```text
observed_i = q + noise_i
E[max(observed_i)] - q ≈ s * sqrt(2 ln N)
```

The default simulation gives:

| N | Mean upward bias | Monte Carlo SEM | Leading approximation |
| ---: | ---: | ---: | ---: |
| 10 | 0.07779 | 0.00175 | 0.10730 |
| 100 | 0.12344 | 0.00123 | 0.15174 |
| 1,000 | 0.15932 | 0.00097 | 0.18585 |
| 10,000 | 0.19343 | 0.00090 | 0.21460 |

These are synthetic outputs, rounded from 300 independent trials per N.
The approximation is a leading term, and also an upper bound on the expected
maximum for centered Gaussian noise. It is loose at small N. Do not present it
as an exact bias correction. `s` is the SD of a board score, not the SD of a
single unit. Monte Carlo SEM in the table measures simulation precision.

Equal true quality makes the failure clear: a larger search wins more noise,
not more quality. Nearly equal quality mixes real gains with this noise.
Correlated candidates reduce the independent-search effect. Heavy tails,
adaptive searches and bounded scores change the approximation. Count all
looks even if the number of independent comparisons is uncertain.

Background: Andrew Nobel's [Gaussian Extreme Values notes](https://nobel.web.unc.edu/wp-content/uploads/sites/13591/2024/04/Gaussian_Extremes.pdf)
give the Gaussian maximum bound and limiting theory. Selection can overfit
an unbiased performance estimate; see Cawley and Talbot (2010),
[On Over-fitting in Model Selection and Subsequent Selection Bias in Performance Evaluation](https://jmlr.org/papers/v11/cawley10a.html), JMLR 11:2079–2107.

## 2. Easier validation than test

Give five candidates a quirk sensitivity `a` from zero to one:

```text
true quality = 0.60 - 0.08*a
validation   = true quality + 0.25*a + small independent noise
test         = true quality          + small independent noise
```

The extra validation term rewards a board-specific shortcut. The strongest
shortcut candidate ranks first on validation and last on test. Its synthetic
scores are about `0.770` and `0.518`, a drop of `0.251`. The printed table shows
all five candidates and their drops. More validation seeds reduce noise;
they do not remove this systematic mismatch.

The example deliberately makes quirk sensitivity trade off against transferable
quality. A constant bias shared by all candidates would shift scores without
reversing their order. This toy mechanism does not diagnose any actual board.

## 3. Sealed selection

Explore 100 candidates on a noisy development board. Freeze three finalists
before opening a separate sealed proxy. Score each finalist once on the same
16 proxy scenarios. Shared scenario effects and independent candidate residuals
give each candidate marginal SD `0.04`. Development and sealed random streams
are independent. Sealing means the finalist identities, code, panel, floor,
look count and acceptance rule are fixed before these outcomes are seen.

For each finalist, compute the sample mean, sample SD and
`SEM = SD / sqrt(16)`. Accept if

```text
mean - SD * sqrt(2 ln k) - 2*SEM > floor
```

Here `k = 3` is the total number of finalists examined on the sealed panel.
The floor is a predeclared, known synthetic value `0.50`. The script compares
this rule with `mean > floor` and `mean - 2*SEM > floor`. It reports the
probability that at least one finalist is accepted in an entire experiment,
not a per-candidate false-positive rate.

Repeat this experiment with all candidate true qualities equal to the floor,
then with gains of `0.075` and `0.12`. In the null experiment every acceptance
is a false win. In the positive experiments a rejection of every finalist is a
missed true win. The stricter rule sharply reduces false wins in this toy
setting. It also misses some moderate gains. Larger gains survive more often.
The output includes counts and rates so the cost is visible.

Default synthetic rates across 1,000 experiments per row:

| Experiment and rate | Mean only | Mean minus 2 SEM | Full lower-bound rule |
| --- | ---: | ---: | ---: |
| No true gain: false win | 0.820 | 0.074 | 0.000 |
| True gain 0.075: missed win | 0.000 | 0.000 | 0.203 |
| True gain 0.120: missed win | 0.000 | 0.000 | 0.000 |

This is a conservative teaching rule, not a calibrated confidence interval or
a distribution-free error guarantee. The requested penalty uses outcome SD,
while uncertainty in a mean scales with SEM. Its SD penalty therefore remains
as the panel grows. No theorem turns this plug-in rule into a chosen error
level under unknown tails or dependence. Zero false wins in 1,000 trials still
allows roughly `3 / 1000` as a one-sided 95% upper rate under independent
Bernoulli experiments. That statement concerns Monte Carlo precision only.

If the floor is an estimated incumbent, compare paired differences with zero
rather than treating its estimate as known. If finalists or rules change after
opening the panel, it is no longer sealed. Log the look and use a fresh panel
for subsequent confirmation. A sealed proxy still cannot fix the systematic
distribution shift from section 2. The repeated trials here each create a
fresh synthetic experiment; they are not repeated chances to pass one panel.

## 4. Paired seeds

Use the same random stream `z_j` for both candidates on each scenario:

```text
A_j = 0.600 + 0.100*z_j
B_j = 0.608 + 0.110*z_j
D_j = B_j - A_j = 0.008 + 0.010*z_j
```

The seed fixes the complete stochastic setup for both candidates. Across 64
synthetic streams, B has the higher mean. Ranking each candidate by its own
`mean - SD` instead picks A. That marginal shrinkage penalizes B's slightly
greater sensitivity to a common nuisance. It answers a different question
from whether B improves on A under matched conditions.

Use the differences `D_j`. Their mean is about `0.00770` and their SEM about
`0.00157`; `mean(D) - 2*SEM(D)` is positive. Treating the two means as independent
gives an SEM about `0.02327` and loses the common-noise cancellation:

```text
Var(B-A) = Var(B) + Var(A) - 2*Cov(A,B)
```

Pairing does not change arithmetic mean order. The inversion here is caused
by separately subtracting marginal SDs. A two-SEM comparison is illustrative,
and still needs selection accounting when many pairs are inspected. Repeated
scorer seeds probe stochastic scoring noise. They do not create independent
biological specimens or establish transfer to a new population.

## 5. Best-of-k boards

Assume a board retains the best valid score among k submissions, one slot
remains, and an additional attempt has no cost or effect on the incumbent.
Under this assumed mechanism, the last valid attempt cannot lower that best
score. This is a decision example, not a statement of current challenge rules.

For a known incumbent score `c` and a possible challenger score `X`, use expected
improvement, also called expected gain here:

```text
EI = EG = E[(X-c)+], where (u)+ = max(u, 0)
E[max(c, X)] = c + EG
```

Estimate that expectation on the same proxy panel. For scenario j, score both
the incumbent `C_j` and challenger `X_j`, then average `max(X_j-C_j, 0)`.
This extends the fixed-c formula to an uncertain incumbent and preserves their
joint variation. Replacing `C_j` with its mean loses dependence. An incumbent
public board number on a different scale is never a valid proxy threshold.
Even matching numeric ranges does not establish matching calibration.

The script builds 256 common proxy scenarios:

```text
C_j = 0.70 + 0.08*z_j
X_j = C_j - 0.01 + 0.04*e_j
near_copy_j = C_j - 0.00001 + 0.00001*e_j
```

The challenger has a lower mean. The lower-bound rule applied to `X-C`, with
floor zero and two compared policies, says keep the incumbent. EG is positive,
so the final free attempt rule says submit the challenger. Its population EG
under this toy Gaussian model is `0.011454`. The script also prints the finite
panel estimate, which differs because of sampling noise.

The near-copy has a high absolute mean but EG near zero. Its changes are tiny
and strongly correlated with the incumbent. An exact copy has EG exactly zero.
The script prints enough decimals to distinguish a small gain from zero.
With zero cost, even a tiny positive EG favors an attempt mathematically; a
predeclared cost or minimum useful gain can justify keeping the incumbent.

For the final slot, choose the eligible challenger with the highest estimated
EG, including keep as zero. Lower bounds address evidence of reliable
improvement. EG addresses expected improvement of a protected maximum. They
optimize different objectives. More remaining slots require considering joint
outcomes across submissions. EG estimates can themselves be overfit if searched
across many candidates. Freeze finalists and the proxy panel first. A positive
proxy EG is not a guarantee of gain on an unknown test board.

## 6. Checklist

- Freeze the decision rule before scoring. Specify the objective, direction,
  floor, finalist limit, uncertainty penalty and tie rule.
- Count every look. Include rejected variants, seed changes, metric choices,
  repeated inspections and changes made after reading a result.
- Keep the panel the same across comparisons. Align scenario IDs, seeds,
  preprocessing and scale. Score the incumbent on that panel too.
- Record what you looked at. Keep candidate identifiers, fixed configurations,
  panel version, seeds, all observed results, rule version and the decision.
  Record discarded candidates as well as the winner.

## Sources and contribution overlap

Rechecked the live [community register](https://api.virtualembryo.ai/challenge/community/register)
on 2026-11-04, including both lists: 15 reviewed contributions and 674 pending
entries. The premise that candidate selection and the winner's curse were
absent did not hold. A listing is not an endorsement or an award decision.
Pending records provide titles and links, without detailed summaries.

| Existing resource | Overlap |
| --- | --- |
| [VEC Winner's Curse Selection Simulator](https://random-guy-05.github.io/vec-participation-workbench/simulations/winner-curse/), pending, filed 2026-10-15 | Direct topic overlap with section 1. No claim to introduce this lesson to the community. |
| [VEC BatchRank](https://github.com/random-guy-05/vec-batchrank), pending | Candidate selection, seed sensitivity and marginal `mean - lambda*SD` ranking. Its public README was inspected. |
| [Reliability Lab 14](https://github.com/random-guy-05/vec-reliability-lab/blob/main/resources/en/14-proxy-rank-reversal.md), pending | Register title explicitly covers reversal between local and official method rankings. |
| [Reliability Lab 16](https://github.com/random-guy-05/vec-reliability-lab/blob/main/resources/en/16-method-family-dependence.md), pending | Register title covers dependence among hyperparameter variants. |
| [VecBench](https://github.com/AJK0921/vecbench), reviewed | Register summary covers scoring noise and holdout ranking reliability. No real outcome numbers are reused. |
| [Local Scoring Wrapper](https://github.com/xxx12e/vec-community-kit/tree/main/vec_local_score), reviewed | Multi-seed local summaries and proxy limitations. |
| [Experiment Plan Card](https://random-guy-05.github.io/vec-participation-workbench/tools/experiment-card/), [RunLog](https://github.com/random-guy-05/vec-runlog) and [Attempt Ledger](https://github.com/random-guy-05/vec-p3-attempt-ledger), pending | Planning, provenance and counting attempts overlap with the checklist. |
| [Paired Median-Difference Bootstrap](https://random-guy-05.github.io/vec-participation-workbench/tools/paired-median-bootstrap/), pending | Register title covers resampling complete pairs; paired inference is established work. |

This independent draft combines synthetic Gaussian maxima, explicit board bias,
sealed acceptance with missed-win costs, common-stream differences and the
protected-maximum EG objective in one executable lesson with tests. Those
components use standard statistical ideas. This bounded review does not
establish global novelty. No existing code, real scientific data or real
outcome tables were copied. The general sources are credited in section 1.
