#!/usr/bin/env python3
"""Synthetic proxy-design experiments.

Standard library only. matplotlib is optional and used only with --plot.
No input files, no network, no challenge data. Every number is simulated.
"""

from __future__ import annotations

import argparse
import math
import random
from dataclasses import dataclass
from statistics import mean, stdev
from typing import Callable, Mapping, Sequence

Vector = list[float]
Observed = Mapping[int, Vector]


# ---------------------------------------------------------------------------
# Small helpers


def mse(pred: Sequence[float], target: Sequence[float]) -> float:
    if not pred or len(pred) != len(target):
        raise ValueError("prediction and target must be nonempty and aligned")
    return sum((p - t) ** 2 for p, t in zip(pred, target)) / len(pred)


def pearson(xs: Sequence[float], ys: Sequence[float]) -> float:
    if len(xs) < 3 or len(xs) != len(ys):
        raise ValueError("need at least three aligned values")
    mx, my = mean(xs), mean(ys)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    if sxx == 0 or syy == 0:
        raise ValueError("correlation is undefined for a constant series")
    return sxy / math.sqrt(sxx * syy)


def slope(xs: Sequence[float], ys: Sequence[float]) -> float:
    """Least-squares slope of ys on xs."""
    mx, my = mean(xs), mean(ys)
    sxx = sum((x - mx) ** 2 for x in xs)
    if sxx == 0:
        raise ValueError("slope is undefined for a constant predictor")
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx


def sem(values: Sequence[float]) -> float:
    return stdev(values) / math.sqrt(len(values)) if len(values) > 1 else float("nan")


# ---------------------------------------------------------------------------
# 1. Horizon leakage


@dataclass(frozen=True)
class Proxy:
    """A local board: predict `target` starting from `anchor`.

    `observed` lists the stage times a method may read. For an extrapolation
    test, nothing after the anchor is observed.
    """

    name: str
    anchor: int
    target: int
    observed: tuple[int, ...]

    @property
    def horizon(self) -> int:
        return self.target - self.anchor


def steps_after_anchor(proxy: Proxy) -> tuple[int, ...]:
    """Observed stage times later than the anchor. Nonempty means leakage."""
    return tuple(t for t in sorted(proxy.observed) if t > proxy.anchor)


def clean_proxy(released: Sequence[int], horizon: int, name: str = "clean") -> Proxy:
    """Latest released stage as target, anchor `horizon` earlier, nothing in between.

    Raises if the released stages cannot support the horizon with at least
    two observed stages at or before the anchor.
    """
    if horizon < 1:
        raise ValueError("horizon must be a positive number of steps")
    times = sorted(set(released))
    target = times[-1]
    anchor = target - horizon
    if anchor not in times:
        raise ValueError(f"no released stage at anchor time {anchor}")
    observed = tuple(t for t in times if t <= anchor)
    if len(observed) < 2:
        raise ValueError("need two observed stages at or before the anchor")
    return Proxy(name, anchor, target, observed)


def simulate_trajectories(
    rng: random.Random,
    genes: int,
    times: Sequence[int],
    slope_sd: float = 1.0,
    curve_sd: float = 0.05,
    noise_sd: float = 0.3,
) -> dict[int, Vector]:
    """Smooth per-gene quadratics, one noisy specimen per stage."""
    params = [
        (rng.gauss(0, 1), rng.gauss(0, slope_sd), rng.gauss(0, curve_sd))
        for _ in range(genes)
    ]
    return {
        t: [a + b * t + c * t * t + rng.gauss(0, noise_sd) for a, b, c in params]
        for t in times
    }


def predict_trend(obs: Observed, proxy: Proxy) -> Vector:
    """Linear extrapolation from the anchor and the stage before it."""
    before = [t for t in proxy.observed if t < proxy.anchor]
    if not before:
        raise ValueError("trend needs a stage before the anchor")
    p, a, h = max(before), proxy.anchor, proxy.horizon
    return [xa + h * (xa - xp) / (a - p) for xa, xp in zip(obs[a], obs[p])]


def predict_reuse(obs: Observed, proxy: Proxy, weight: float = 0.5) -> Vector:
    """Trend blended with the observed stage nearest to the target.

    This is the kind of shortcut that copies whatever is closest. It does
    not check whether the nearest stage would exist at test time.
    """
    candidates = [t for t in proxy.observed if t != proxy.target]
    nearest = min(candidates, key=lambda t: (abs(proxy.target - t), -t))
    trend = predict_trend(obs, proxy)
    return [(1 - weight) * x + weight * n for x, n in zip(trend, obs[nearest])]


def horizon_boards() -> list[Proxy]:
    released = (0, 1, 2, 3, 4)
    return [
        Proxy("leaky proxy", anchor=2, target=4, observed=(0, 1, 2, 3)),
        clean_proxy(released, horizon=2, name="clean proxy"),
        Proxy("test", anchor=4, target=6, observed=released),
    ]


def horizon_leakage(repeats: int = 40, genes: int = 400, seed: int = 101) -> list[dict]:
    rng = random.Random(seed)
    boards = horizon_boards()
    diffs: dict[str, list[float]] = {b.name: [] for b in boards}
    trend_err: dict[str, list[float]] = {b.name: [] for b in boards}
    reuse_err: dict[str, list[float]] = {b.name: [] for b in boards}
    for _ in range(repeats):
        obs = simulate_trajectories(rng, genes, range(7))
        for b in boards:
            target = obs[b.target]
            visible = {t: obs[t] for t in b.observed}
            et = mse(predict_trend(visible, b), target)
            er = mse(predict_reuse(visible, b), target)
            trend_err[b.name].append(et)
            reuse_err[b.name].append(er)
            diffs[b.name].append(er - et)
    return [
        {
            "board": b.name,
            "anchor": b.anchor,
            "target": b.target,
            "horizon": b.horizon,
            "after_anchor": steps_after_anchor(b),
            "trend_mse": mean(trend_err[b.name]),
            "reuse_mse": mean(reuse_err[b.name]),
            "reuse_minus_trend": mean(diffs[b.name]),
            "sem": sem(diffs[b.name]),
        }
        for b in boards
    ]


# ---------------------------------------------------------------------------
# 3. Step persistence: a diagnostic, and the coefficient for a last-step rule


def simulate_steps(
    rng: random.Random,
    genes: int,
    stages: int,
    persistence: float,
    noise_sd: float = 0.2,
) -> dict[int, Vector]:
    """Stage-to-stage change = shared per-gene drift + a new program each step.

    Each true step has unit variance per gene. `persistence` is the share of
    that variance carried by the drift, so it is the true correlation between
    any two different steps, at every lag. Observations add one-specimen noise.
    """
    if not 0 <= persistence <= 1:
        raise ValueError("persistence must be between 0 and 1")
    keep, new = math.sqrt(persistence), math.sqrt(1 - persistence)
    drift = [rng.gauss(0, 1) for _ in range(genes)]
    truth = [rng.gauss(0, 1) for _ in range(genes)]
    obs = {0: [x + rng.gauss(0, noise_sd) for x in truth]}
    for k in range(1, stages):
        truth = [x + keep * d + new * rng.gauss(0, 1) for x, d in zip(truth, drift)]
        obs[k] = [x + rng.gauss(0, noise_sd) for x in truth]
    return obs


def simulate_rotating(
    rng: random.Random,
    genes: int,
    stages: int,
    noise_sd: float = 0.2,
) -> dict[int, Vector]:
    """Counterexample where dependence is not equal across lags.

    The true step from stage k to k+1 is A*cos(k*pi/4) + B*sin(k*pi/4), with
    independent standard-normal A and B per gene. Steps one apart correlate at
    cos(pi/4) = 0.71; steps two apart at cos(pi/2) = 0.
    """
    coef = [(rng.gauss(0, 1), rng.gauss(0, 1)) for _ in range(genes)]
    truth = [rng.gauss(0, 1) for _ in range(genes)]
    obs = {0: [x + rng.gauss(0, noise_sd) for x in truth]}
    for k in range(1, stages):
        c, s = math.cos((k - 1) * math.pi / 4), math.sin((k - 1) * math.pi / 4)
        truth = [x + a * c + b * s for x, (a, b) in zip(truth, coef)]
        obs[k] = [x + rng.gauss(0, noise_sd) for x in truth]
    return obs


def step(obs: Observed, start: int, end: int) -> Vector:
    return [b - a for a, b in zip(obs[start], obs[end])]


def transfer_check(obs: Observed, last: int) -> dict:
    """Step correlations on released stages up to `last`.

    `disjoint_r` correlates (last-3 -> last-2) with (last-1 -> last). The two
    steps share no stage, so it measures lag-two persistence of the true
    steps; independent specimen noise adds no covariance but attenuates the
    observed correlation. `adjacent_r` correlates (last-2 -> last-1) with (last-1 -> last).
    `shrink` is the adjacent least-squares slope, clipped at zero. A last-step
    rule predicts one step from the step right before it, so its coefficient
    must come from an adjacent pair: the shared stage's noise is part of the
    rule's own error.
    """
    early = step(obs, last - 3, last - 2)
    late = step(obs, last - 1, last)
    adjacent_a = step(obs, last - 2, last - 1)
    return {
        "disjoint_r": pearson(early, late),
        "adjacent_r": pearson(adjacent_a, late),
        "shrink": max(0.0, slope(adjacent_a, late)),
    }


def rolling_origin_gain(obs: Observed, last: int) -> float:
    """Copy MSE minus shrunk-rule MSE, one step ahead, on released stages only.

    For each origin o with two earlier stages and o + 1 <= last, fit the
    coefficient on the adjacent pair (o-2 -> o-1, o-1 -> o), predict o + 1
    from stages up to o, and score against the released stage o + 1.
    Positive means the rule beat copying on the stages you have.
    """
    gains = []
    for o in range(2, last):
        beta = max(0.0, slope(step(obs, o - 2, o - 1), step(obs, o - 1, o)))
        last_step = step(obs, o - 1, o)
        pred = [a + beta * s for a, s in zip(obs[o], last_step)]
        gains.append(mse(obs[o], obs[o + 1]) - mse(pred, obs[o + 1]))
    if not gains:
        raise ValueError("need at least four released stages")
    return mean(gains)


@dataclass(frozen=True)
class StepCase:
    name: str
    persistence: float | None  # None: rotating counterexample
    noise_sd: float


STEP_CASES = (
    StepCase("p=0.80 s=0.2", 0.8, 0.2),
    StepCase("p=0.05 s=0.2", 0.05, 0.2),
    StepCase("p=0.30 s=0.6", 0.3, 0.6),
    StepCase("rotating s=0.2", None, 0.2),
)


def disjoint_transfer(
    cases: Sequence[StepCase] = STEP_CASES,
    genes: int = 800,
    repeats: int = 20,
    seed: int = 202,
) -> list[dict]:
    """Released stages 0..4, test stage 5. Means over repeats; sd is across repeats."""
    rng = random.Random(seed)
    rows = []
    for case in cases:
        draws = []
        for _ in range(repeats):
            if case.persistence is None:
                obs = simulate_rotating(rng, genes, stages=6, noise_sd=case.noise_sd)
            else:
                obs = simulate_steps(rng, genes, 6, case.persistence, case.noise_sd)
            released_last, test = 4, 5
            check = transfer_check(obs, released_last)
            last_step = step(obs, released_last - 1, released_last)
            anchor, target = obs[released_last], obs[test]
            beta = check["shrink"]
            copy = mse(anchor, target)
            shrunk = mse([a + beta * s for a, s in zip(anchor, last_step)], target)
            draws.append(
                {
                    **check,
                    "beta": beta,
                    "copy_mse": copy,
                    "extrapolate_mse": mse([a + s for a, s in zip(anchor, last_step)], target),
                    "shrunk_mse": shrunk,
                    "gain": copy - shrunk,
                    "released_gain": rolling_origin_gain(obs, released_last),
                }
            )
        row = {"case": case.name, "persistence": case.persistence, "noise_sd": case.noise_sd}
        for key in draws[0]:
            row[key] = mean(d[key] for d in draws)
        row["gain_sd"] = stdev(d["gain"] for d in draws)
        rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# 4. When two proxies disagree


@dataclass(frozen=True)
class Design:
    """What a board looks like, known before any candidate is scored."""

    name: str
    horizon: int
    depth_step: float  # log depth change in the last observed step
    specimens: int     # specimens averaged per stage


def criteria_met(test: Design, proxy: Design, depth_tol: float = 0.15) -> tuple[bool, bool, bool]:
    """Horizon equal, depth change within `depth_tol`, specimen count equal.

    The tolerance is an illustrative choice. Declare it before scoring.
    """
    return (
        proxy.horizon == test.horizon,
        abs(proxy.depth_step - test.depth_step) <= depth_tol,
        proxy.specimens == test.specimens,
    )


def design_screen(test: Design, proxies: Sequence[Design], depth_tol: float = 0.15) -> Design | None:
    """Conservative design screen. Takes no scores.

    Returns a proxy only if exactly one proxy meets all three criteria.
    Otherwise returns None: no proxy settles the question, so build a
    matched one or report that the question is open.
    """
    full = [d for d in proxies if all(criteria_met(test, d, depth_tol))]
    return full[0] if len(full) == 1 else None


def count_matches_choice(test: Design, proxies: Sequence[Design], depth_tol: float = 0.15) -> Design | None:
    """The rule this tutorial does not recommend: most criteria met wins.

    Kept to show its failure. Ties and zero matches return None.
    """
    counts = [(sum(criteria_met(test, d, depth_tol)), d) for d in proxies]
    counts.sort(key=lambda c: c[0], reverse=True)
    if not counts or counts[0][0] == 0:
        return None
    if len(counts) > 1 and counts[0][0] == counts[1][0]:
        return None
    return counts[0][1]


def simulate_board(
    rng: random.Random,
    design: Design,
    genes: int,
    persistence: float = 0.65,
    noise_sd: float = 0.3,
) -> float:
    """Gain of 'extrapolate the last step' over 'copy the anchor'. Positive = helps.

    Depth is a stand-in: a stage sequenced shallower is shifted by
    -sensitivity_g * |log depth change|. This is not a count model.
    """
    keep, new = math.sqrt(persistence), math.sqrt(1 - persistence)
    noise = noise_sd / math.sqrt(design.specimens)
    h = design.horizon
    sens = [abs(rng.gauss(0.5, 0.25)) for _ in range(genes)]
    drift = [rng.gauss(0, 1) for _ in range(genes)]
    x_prev = [rng.gauss(0, 1) for _ in range(genes)]
    x_anchor = [x + keep * d + new * rng.gauss(0, 1) for x, d in zip(x_prev, drift)]
    x_target = list(x_anchor)
    for _ in range(h):
        x_target = [x + keep * d + new * rng.gauss(0, 1) for x, d in zip(x_target, drift)]
    # The anchor stage is shallower than the stage before it; the target is
    # scored at the anchor's depth, so only the last observed step carries depth.
    shift = [design.depth_step * s for s in sens]
    o_prev = [x + rng.gauss(0, noise) for x in x_prev]
    o_anchor = [x + s + rng.gauss(0, noise) for x, s in zip(x_anchor, shift)]
    o_target = [x + s + rng.gauss(0, noise) for x, s in zip(x_target, shift)]
    copy = o_anchor
    extrap = [a + h * (a - p) for a, p in zip(o_anchor, o_prev)]
    return mse(copy, o_target) - mse(extrap, o_target)


def proxy_disagreement(repeats: int = 40, genes: int = 300, seed: int = 303) -> dict:
    test = Design("test", horizon=2, depth_step=-0.8, specimens=1)
    proxies = [
        Design("proxy A", horizon=1, depth_step=0.0, specimens=3),
        Design("proxy B", horizon=2, depth_step=-0.7, specimens=1),
    ]
    # Freeze the choice before scoring anything.
    chosen = design_screen(test, proxies)
    rng = random.Random(seed)
    gains = {d.name: [simulate_board(rng, d, genes) for _ in range(repeats)] for d in (*proxies, test)}
    test_sign = math.copysign(1, mean(gains["test"]))
    rows = []
    for d in (*proxies, test):
        g = gains[d.name]
        rows.append(
            {
                "board": d.name,
                "horizon": d.horizon,
                "depth_step": d.depth_step,
                "specimens": d.specimens,
                "gain": mean(g),
                "sd": stdev(g),
                "sign_matches_test": sum(math.copysign(1, x) == test_sign for x in g) / len(g),
            }
        )
    # One factor at a time, starting from proxy A and moving toward the test.
    a = proxies[0]
    variants = [
        Design("A, test horizon", test.horizon, a.depth_step, a.specimens),
        Design("A, test depth", a.horizon, test.depth_step, a.specimens),
        Design("A, test specimens", a.horizon, a.depth_step, test.specimens),
    ]
    ablation = []
    for d in variants:
        g = [simulate_board(rng, d, genes) for _ in range(repeats)]
        ablation.append({"board": d.name, "gain": mean(g), "sem": sem(g)})
    # Counterexample: counting matches picks a proxy with the wrong sign.
    c = Design("proxy C", horizon=2, depth_step=0.0, specimens=1)
    d = Design("proxy D", horizon=1, depth_step=-0.8, specimens=2)
    counter = []
    for x in (c, d):
        g = [simulate_board(rng, x, genes) for _ in range(repeats)]
        counter.append({"board": x.name, "horizon": x.horizon, "depth_step": x.depth_step,
                        "specimens": x.specimens, "matches": sum(criteria_met(test, x)),
                        "gain": mean(g), "sem": sem(g)})
    count_pick = count_matches_choice(test, [c, d])
    screen_pick = design_screen(test, [c, d])
    return {
        "chosen": chosen.name if chosen else None,
        "rows": rows,
        "ablation": ablation,
        "counter": counter,
        "counter_count_pick": count_pick.name if count_pick else None,
        "counter_screen_pick": screen_pick.name if screen_pick else None,
    }


# ---------------------------------------------------------------------------
# Output


def report(repeats: int, genes: int) -> str:
    out = []
    out.append("1. Horizon leakage (mean squared error, lower is better)")
    out.append(f"{'board':<12} {'anchor->target':>14} {'after anchor':>13} {'trend':>8} {'reuse':>8} {'reuse-trend':>12} {'sem':>6}")
    for r in horizon_leakage(repeats=repeats, genes=genes):
        span = f"{r['anchor']}->{r['target']}"
        after = ",".join(map(str, r["after_anchor"])) or "none"
        out.append(
            f"{r['board']:<12} {span:>14} {after:>13} {r['trend_mse']:8.3f} {r['reuse_mse']:8.3f}"
            f" {r['reuse_minus_trend']:12.3f} {r['sem']:6.3f}"
        )
    out.append("")
    out.append("2. Depth mismatch: covered by an earlier filing; see README.")
    out.append("")
    out.append("3. Step persistence (released stages 0..4, test stage 5; means over repeats)")
    out.append(
        f"{'case':<15} {'disjoint r':>10} {'adjacent r':>10} {'beta':>6} {'copy':>6} {'extrap':>6}"
        f" {'shrunk':>6} {'gain':>7} {'sd':>6} {'released':>8}"
    )
    for r in disjoint_transfer(genes=max(genes * 2, 200), repeats=max(repeats // 2, 2)):
        out.append(
            f"{r['case']:<15} {r['disjoint_r']:10.3f} {r['adjacent_r']:10.3f} {r['beta']:6.3f}"
            f" {r['copy_mse']:6.3f} {r['extrapolate_mse']:6.3f} {r['shrunk_mse']:6.3f}"
            f" {r['gain']:7.3f} {r['gain_sd']:6.3f} {r['released_gain']:8.3f}"
        )
    out.append("")
    res = proxy_disagreement(repeats=repeats, genes=genes)
    out.append("4. Two proxies disagree (gain of extrapolating over copying; positive = helps)")
    out.append(f"{'board':<8} {'horizon':>7} {'depth':>6} {'specimens':>9} {'gain':>7} {'sd':>6} {'sign as test':>12}")
    for r in res["rows"]:
        out.append(
            f"{r['board']:<8} {r['horizon']:7d} {r['depth_step']:6.2f} {r['specimens']:9d}"
            f" {r['gain']:7.3f} {r['sd']:6.3f} {r['sign_matches_test']:12.2f}"
        )
    out.append(f"Design screen, before scoring: {res['chosen'] or 'abstain'}")
    out.append("One factor changed from proxy A toward the test:")
    for r in res["ablation"]:
        out.append(f"  {r['board']:<18} gain {r['gain']:7.3f}  sem {r['sem']:.3f}")
    out.append("Counterexample for counting matches (same test):")
    for r in res["counter"]:
        out.append(
            f"  {r['board']:<8} horizon {r['horizon']} depth {r['depth_step']:5.2f} specimens {r['specimens']}"
            f"  matches {r['matches']}  gain {r['gain']:7.3f}  sem {r['sem']:.3f}"
        )
    out.append(f"  counting matches picks: {res['counter_count_pick'] or 'none'}")
    out.append(f"  design screen picks:    {res['counter_screen_pick'] or 'abstain'}")
    return "\n".join(out)


def plot(path: str, repeats: int, genes: int) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise SystemExit("matplotlib is not installed; run without --plot") from exc
    rows = horizon_leakage(repeats=repeats, genes=genes)
    rng = random.Random(202)
    obs = simulate_steps(rng, 1000, stages=6, persistence=0.8)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3.6))
    ax1.bar([r["board"] for r in rows], [r["reuse_minus_trend"] for r in rows], color="#4a6fa5")
    ax1.axhline(0, color="black", linewidth=0.8)
    ax1.set_ylabel("reuse MSE minus trend MSE")
    ax1.set_title("Below zero: reuse looks better")
    ax2.scatter(step(obs, 1, 2), step(obs, 3, 4), s=3, alpha=0.4, color="#4a6fa5")
    ax2.set_xlabel("step 1->2")
    ax2.set_ylabel("step 3->4")
    ax2.set_title("Disjoint steps, persistence 0.8")
    fig.tight_layout()
    fig.savefig(path, dpi=120)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repeats", type=int, default=40, help="synthetic repeats per board")
    parser.add_argument("--genes", type=int, default=400, help="synthetic features per stage")
    parser.add_argument("--plot", metavar="PNG", help="also save a figure (needs matplotlib)")
    args = parser.parse_args(argv)
    if args.repeats < 2 or args.genes < 10:
        parser.error("use --repeats >= 2 and --genes >= 10")
    print(report(args.repeats, args.genes))
    if args.plot:
        plot(args.plot, args.repeats, args.genes)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
