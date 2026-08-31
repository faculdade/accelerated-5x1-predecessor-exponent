#!/usr/bin/env python3
"""
E-002, part 8: does the reading differ between basin roots and the rest?

The paper says Conjecture 7.2 "carries almost all of the weight", because only
18 of the 300 sampled roots reach a known cycle and the other 282 are disjoint
from the set Volkov's Q-tilde counts. That is the paper's account of its weakest
link, and round 19 (R19-05) found it stops one step short: the split it worries
about is measurable on the committed data, and nobody measured it.

The 18 basin roots come from `cycle_membership.py` (recomputed here, not copied),
their per-root readings are in `data/q5_arith_b15.txt`, and splitting the sample
costs nothing. If the two subsamples read the same, that is direct empirical
support for exactly the transfer Conjecture 7.2 is invoked to make.

Estimator: `summary.py`'s own per-root slope and Aitken, imported. The pooled
value over all 300 roots must reproduce `out/summary_b15_d10.log`'s 0.64926;
the script checks that and aborts if not.

Serial (Rule 9c): resampling precomputed per-root curves, about 1.3 s.

Run: python3 basin_split.py
"""
import random
import sys

sys.path.insert(0, ".")
from summary import load, slope, aitken, mean            # noqa: E402
from cycle_membership import sample_e002, CYCLE5, STEPS  # noqa: E402

DEC, SEED, B = 10, 17, 4000
POOLED = 0.64926


def basin_roots():
    """the sampled roots whose forward orbit reaches a known cycle"""
    out = set()
    for u in sample_e002():
        v = u
        for _ in range(STEPS):
            if v in CYCLE5:
                out.add(u)
                break
            v = 5 * v + 1
            while v % 2 == 0:
                v //= 2
    return out


def main():
    _, roots, mats, cp_lo, buf_lo, n_cp, n_buf, _ = load("data/q5_arith_b15.txt")
    ci = DEC - 1 - cp_lo
    use = [bi for bi in range(n_buf) if buf_lo + bi > DEC]
    curves = {}
    for root, m in zip(roots, mats):
        r = [slope(m, ci, ci + 1, bi) for bi in use]
        if all(v is not None for v in r):
            curves[root] = r

    def est(rs):
        rows = [curves[r] for r in rs if r in curves]
        if len(rows) < 3:
            return None
        return aitken([mean([x[b] for x in rows]) for b in range(len(use))])

    allr = sorted(curves)
    pooled = est(allr)
    if round(pooled, 5) != POOLED:
        raise SystemExit(f"ERROR: pooled reads {pooled:.5f}, "
                         f"out/summary_b15_d10.log says {POOLED:.5f}")

    basin = sorted(basin_roots() & set(curves))
    esc = sorted(set(allr) - set(basin))
    b_est, e_est = est(basin), est(esc)
    print(f"arithmetic tree, decade 1e{DEC-1} -> 1e{DEC}, b15 grid\n")
    print(f"  all {len(allr)} roots           {pooled:.5f}   "
          f"(reproduces out/summary_b15_d10.log)")
    print(f"  the {len(basin)} basin roots      {b_est:.5f}   "
          f"these are the roots Volkov's Q-tilde counts")
    print(f"  the {len(esc)} remaining roots {e_est:.5f}   (not proven to escape; none is)")
    print(f"  difference             {b_est - e_est:+.5f}")

    rng = random.Random(SEED)
    diffs = []
    for _ in range(B):
        bb = [basin[rng.randrange(len(basin))] for _ in basin]
        ee = [esc[rng.randrange(len(esc))] for _ in esc]
        x, y = est(bb), est(ee)
        if x is not None and y is not None:
            diffs.append(x - y)
    diffs.sort()
    lo, hi = diffs[int(0.025 * len(diffs))], diffs[int(0.975 * len(diffs)) - 1]
    print(f"  95% bootstrap interval [{lo:+.5f}, {hi:+.5f}]")

    # the basin arm's own interval, which is what says how much this can settle
    bb = []
    for _ in range(B):
        s = [basin[rng.randrange(len(basin))] for _ in basin]
        v = est(s)
        if v is not None:
            bb.append(v)
    bb.sort()
    blo, bhi = bb[int(0.025 * len(bb))], bb[int(0.975 * len(bb)) - 1]
    print(f"  the basin arm alone: {b_est:.5f}, 95% [{blo:.5f}, {bhi:.5f}], "
          f"half-width {(bhi-blo)/2:.5f}")

    # R20-06: an earlier version of this block printed the distance from
    # eta_5,BP as though it meant something, and called failure to reject
    # "direct support". Both are wrong and the numbers here say why.
    print(f"\n  WHAT THIS DOES AND DOES NOT SHOW. The test does not separate the")
    print(f"  two subsamples, and at n={len(basin)} it could not: the interval on the")
    print(f"  difference, [{lo:+.5f}, {hi:+.5f}], admits values larger than the")
    print(f"  calibration band itself (0.00371). Failure to separate is weaker")
    print(f"  than agreement. The basin arm's own half-width, {(bhi-blo)/2:.5f}, is")
    print(f"  larger than its distance from eta_5,BP, so that distance carries no")
    print(f"  information either. What can be said: nothing in this sample")
    print(f"  contradicts the transfer Conjecture 7.2 is invoked to make.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
