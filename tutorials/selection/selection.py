#!/usr/bin/env python3
"""Synthetic selection experiments. Standard library only; no input files or network."""

from __future__ import annotations

import argparse
import math
import random
from dataclasses import dataclass
from statistics import mean, stdev
from typing import Sequence


@dataclass(frozen=True)
class Summary:
    mean: float
    sd: float
    sem: float


def summarize(values: Sequence[float]) -> Summary:
    if len(values) < 2 or not all(math.isfinite(x) for x in values):
        raise ValueError("at least two finite panel outcomes are required")
    sd = stdev(values)
    return Summary(mean(values), sd, sd / math.sqrt(len(values)))


def extreme_scale(sd: float, count: int) -> float:
    if not isinstance(count, int) or count < 1:
        raise ValueError("count must be a positive integer")
    if not math.isfinite(sd) or sd < 0:
        raise ValueError("sd must be finite and nonnegative")
    return sd * math.sqrt(2 * math.log(count))


def lower_bound(values: Sequence[float], k: int) -> float:
    """Conservative teaching rule, not a calibrated confidence bound."""
    s = summarize(values)
    return s.mean - extreme_scale(s.sd, k) - 2 * s.sem


def paired_differences(
    challenger: Sequence[float], incumbent: Sequence[float]
) -> list[float]:
    if not challenger or len(challenger) != len(incumbent):
        raise ValueError("use nonempty outcomes from the same aligned panel")
    if not all(math.isfinite(x) for x in (*challenger, *incumbent)):
        raise ValueError("panel outcomes must be finite")
    return [x - c for x, c in zip(challenger, incumbent)]


def expected_gain(
    challenger: Sequence[float], incumbent: Sequence[float]
) -> float:
    """Empirical E[(X-C)+] over matched proxy scenarios; also called EI."""
    return mean(max(d, 0.0) for d in paired_differences(challenger, incumbent))


def winners_curse(repeats: int = 300, seed: int = 11) -> list[dict]:
    """Each draw is the exact Gaussian distribution of a small board mean."""
    if repeats < 2:
        raise ValueError("repeats must be at least two")
    rng = random.Random(seed)
    # A board has 16 iid units, each with noise SD 0.20.
    board_sd = 0.20 / math.sqrt(16)
    quality = 0.50
    rows = []
    for n in (10, 100, 1_000, 10_000):
        optimism = []
        for _ in range(repeats):
            best = max(quality + rng.gauss(0, board_sd) for _ in range(n))
            optimism.append(best - quality)
        s = summarize(optimism)
        rows.append(dict(n=n, bias=s.mean, mc_sem=s.sem,
                         approximation=extreme_scale(board_sd, n)))
    return rows


def validation_quirk(seed: int = 13) -> list[dict]:
    rng = random.Random(seed)
    rows = []
    for index in range(5):
        quirk = index / 4
        quality = 0.60 - 0.08 * quirk
        validation = quality + 0.25 * quirk + rng.gauss(0, 0.003)
        test = quality + rng.gauss(0, 0.003)
        rows.append(dict(candidate=f"candidate-{index}", quirk=quirk,
                         quality=quality, validation=validation, test=test,
                         drop=validation - test))
    return rows


def sealed_selection(
    effect: float, trials: int = 1000, seed: int = 19, k: int = 3
) -> dict:
    """Fresh development and sealed boards in each independent toy trial.

    All 100 candidates have quality floor + effect. Exploration chooses k.
    The independent sealed panel is opened once for those finalists only.
    """
    if trials < 1 or not 1 <= k <= 100 or not isinstance(k, int):
        raise ValueError("positive trials and integer k in [1, 100] required")
    if not math.isfinite(effect):
        raise ValueError("effect must be finite")
    development = random.Random(seed)
    sealed = random.Random(seed + 1)
    floor = 0.50
    accepted = dict(mean_only=0, two_sem=0, lower_bound=0)
    for _ in range(trials):
        scores = [floor + effect + development.gauss(0, 0.08)
                  for _ in range(100)]
        finalists = sorted(range(100), key=scores.__getitem__, reverse=True)[:k]
        # The same 16 scenario IDs and common effects serve every finalist.
        # Total marginal outcome variance is 0.04**2 for each candidate.
        common = [sealed.gauss(0, 0.02) for _ in range(16)]
        residual_sd = math.sqrt(0.04**2 - 0.02**2)
        panels = {i: [floor + effect + nuisance + sealed.gauss(0, residual_sd)
                      for nuisance in common] for i in finalists}
        stats = [summarize(panel) for panel in panels.values()]
        accepted["mean_only"] += any(s.mean > floor for s in stats)
        accepted["two_sem"] += any(s.mean - 2 * s.sem > floor for s in stats)
        accepted["lower_bound"] += any(lower_bound(p, k) > floor
                                       for p in panels.values())
    rates = {name: count / trials for name, count in accepted.items()}
    return dict(effect=effect, trials=trials, k=k, accepted=accepted,
                acceptance=rates,
                missed={name: 1 - rate for name, rate in rates.items()})


def paired_seeds(seed: int = 17) -> dict:
    rng = random.Random(seed)
    stream = [rng.gauss(0, 1) for _ in range(64)]
    a = [0.600 + 0.100 * z for z in stream]
    b = [0.608 + 0.110 * z for z in stream]
    sa, sb = summarize(a), summarize(b)
    difference = summarize(paired_differences(b, a))
    return dict(a=sa, b=sb, difference=difference,
                unpaired_sem=math.sqrt(sa.sem**2 + sb.sem**2),
                shrink_a=sa.mean - sa.sd, shrink_b=sb.mean - sb.sd)


def best_of_k(seed: int = 23) -> dict:
    rng = random.Random(seed)
    common = [rng.gauss(0, 1) for _ in range(256)]
    extra = [rng.gauss(0, 1) for _ in common]
    incumbent = [0.70 + 0.08 * z for z in common]
    challenger = [c - 0.01 + 0.04 * e for c, e in zip(incumbent, extra)]
    near_copy = [c - 0.00001 + 0.00001 * e for c, e in zip(incumbent, extra)]
    rows = []
    # k counts the two compared policies here, not a claim about board rules.
    for name, panel in (("challenger", challenger), ("near-copy", near_copy)):
        difference = paired_differences(panel, incumbent)
        rows.append(dict(candidate=name, mean=mean(panel),
                         mean_difference=mean(difference),
                         lower_bound=lower_bound(difference, k=2),
                         eg=expected_gain(panel, incumbent)))
    # Independent Gaussian residual: exact population EG for the challenger.
    delta, sd = -0.01, 0.04
    z = delta / sd
    analytic_eg = sd * math.exp(-z*z / 2) / math.sqrt(2 * math.pi)
    analytic_eg += delta * (1 + math.erf(z / math.sqrt(2))) / 2
    return dict(incumbent_mean=mean(incumbent), rows=rows,
                analytic_eg=analytic_eg)


def positive_int(value: str) -> int:
    number = int(value)
    if number < 2:
        raise argparse.ArgumentTypeError("use an integer of at least two")
    return number


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repeats", type=positive_int, default=300,
                        help="independent winner's-curse boards (default: 300)")
    parser.add_argument("--trials", type=positive_int, default=1000,
                        help="independent sealed experiments (default: 1000)")
    args = parser.parse_args(argv)
    print("SYNTHETIC ONLY. Arbitrary higher-is-better units. Fixed toy seeds.")
    print("\n1. Winner's curse: true quality 0.50; board noise SD 0.05")
    print("N       bias     MC_SEM   sd*sqrt(2 ln N)")
    for row in winners_curse(args.repeats):
        print(f"{row['n']:5d}   {row['bias']:.5f}  {row['mc_sem']:.5f}  "
              f"{row['approximation']:.5f}")
    print("\n2. Easier validation than test")
    print("candidate    quirk  quality validation test    drop")
    for row in validation_quirk():
        print(f"{row['candidate']:12s} {row['quirk']:.2f}   {row['quality']:.3f}   "
              f"{row['validation']:.3f}     {row['test']:.3f}  {row['drop']:.3f}")
    print("\n3. Sealed selection: 100 explored, 3 finalists, 16 sealed units")
    print("effect   rule         accepted/trials   false-win or missed-win rate")
    for effect in (0.0, 0.075, 0.12):
        result = sealed_selection(effect, args.trials)
        for rule, count in result['accepted'].items():
            rate = (result['acceptance'] if effect == 0 else result['missed'])[rule]
            label = "false-win" if effect == 0 else "missed-win"
            print(f"{effect:.3f}    {rule:12s} {count:4d}/{args.trials:<6d}    "
                  f"{label}={rate:.3f}")
    print("Zero observed false wins is not a zero-risk guarantee.")
    print("\n4. Paired seeds: compare B - A on 64 common streams")
    paired = paired_seeds()
    for name in ('a', 'b'):
        s = paired[name]
        print(f"{name.upper()}: mean={s.mean:.5f} SD={s.sd:.5f} "
              f"mean-SD={paired['shrink_' + name]:.5f}")
    d = paired['difference']
    print(f"B-A: mean={d.mean:.5f} paired_SEM={d.sem:.5f} "
          f"unpaired_SEM={paired['unpaired_sem']:.5f} "
          f"mean-2SEM={d.mean - 2*d.sem:.5f}")
    print("\n5. Best-of-k: last free attempt; incumbent on the SAME proxy panel")
    best = best_of_k()
    print(f"incumbent mean={best['incumbent_mean']:.5f}")
    for row in best['rows']:
        decision = "submit" if row['eg'] > 0 else "keep"
        bound_decision = "submit" if row['lower_bound'] > 0 else "keep"
        print(f"{row['candidate']}: mean={row['mean']:.5f} "
              f"mean_delta={row['mean_difference']:.6f} "
              f"lower_bound_delta={row['lower_bound']:.6f} "
              f"EG={row['eg']:.8f}; lower-bound={bound_decision}, EG={decision}")
    print(f"challenger population EG={best['analytic_eg']:.6f}")
    print("\n6. Checklist: freeze rule; count every look; keep panel fixed; log every look.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
