#!/usr/bin/env python3
"""
E-002, part 6: what the calibration band and the separation are actually worth.

Round 18 added a sentence to the paper saying that because the constructions
share one root list, each band is a paired difference in which the between-root
variation cancels. Round 19 (R19-02) measured that and it is false. This script
is the measurement, so the claim and its refutation both live in the deposit.

Pairing cancels variance only when the quantities are CORRELATED across modes.
They are not: at a checkpoint decade of 1e9 to 1e10 a stochastic construction's
tree has forgotten its root, so the spread of its per-root reading is its own
realization noise rather than a common function of the root. The script prints the per-root correlations and, beside EVERY pairwise
difference, the value the standard deviation would take if the two sides were
independent. Most of the 21 pairs agree closely; a few narrow, the largest by
22% between the two constructions that draw an independent residue and share a
realization. The two pairs the paper quotes against the arithmetic tree have
slightly NEGATIVE correlation, so joint resampling widens them rather than
narrowing. "Every pair is unaffected" would be too strong, and rounds 21 and 23
caught the paper saying versions of it.

Consequence for the paper. The wider band is a spread of the point readings
and carries a wide interval of its own, so a ratio to it ("N band widths") is
not a two-significant-figure quantity. The SEPARATION between the arithmetic
tree and the construction retuned to 0.678 is tight, and that is what the paper
now quotes.

The estimator is `summary.py`'s own `decade_estimate`, imported rather than
reimplemented, so the point readings here must equal `out/summary_b15_d10.log`;
the script checks that and aborts if not.

Serial (Rule 9c): 4000 joint resamples over 300 roots and all seven datasets is arithmetic
on precomputed per-root curves and takes about 0.9 s, far under the one-minute
line at which .claude/rules/experiments.md asks for parallel decomposition.

Run: python3 paired_bootstrap.py
"""
import itertools
import math
import random
import statistics
import sys

sys.path.insert(0, ".")
from summary import load, slope, aitken, mean          # noqa: E402

MODES = ["arith", "cyc", "cycq500", "cycq505", "iid", "iidq500", "iidq505"]
# the two constructions built to 0.678; the conservative separation is measured
# to whichever of them reads lower in a given resample (round 22, external)
HIGH = ["cycq505", "iidq505"]
DEC, SEED, B = 10, 17, 4000
# the three constructions built to 0.650919: the integer recursion, its
# real-valued relaxation, and the unconstrained one
WIDE = ["cyc", "cycq500", "iid", "iidq500"]   # all four built to 0.650919
EXPECTED = {"arith": 0.64926, "cyc": 0.65122, "cycq500": 0.64981,
            "cycq505": 0.67748, "iid": 0.64751,
            "iidq500": 0.64738, "iidq505": 0.67358}


def curves(mode):
    _, roots, mats, cp_lo, buf_lo, n_cp, n_buf, _ = load(f"data/q5_{mode}_b15.txt")
    ci = DEC - 1 - cp_lo
    use = [bi for bi in range(n_buf) if buf_lo + bi > DEC]
    out = {}
    for root, m in zip(roots, mats):
        r = [slope(m, ci, ci + 1, bi) for bi in use]
        if all(v is not None for v in r):
            out[root] = r
    return out, len(use)


def main():
    data, nb = {}, None
    for m in MODES:
        data[m], nb = curves(m)
    common = sorted(set.intersection(*(set(v) for v in data.values())))

    def point(mode, idx):
        rows = [data[mode][common[i]] for i in idx]
        return aitken([mean([r[b] for r in rows]) for b in range(nb)])

    base = list(range(len(common)))
    pts = {m: point(m, base) for m in MODES}
    for m in MODES:
        if round(pts[m], 5) != EXPECTED[m]:
            raise SystemExit(f"ERROR: {m} reads {pts[m]:.5f}, out/summary_b15_d10.log "
                             f"says {EXPECTED[m]:.5f}; this is not the paper's estimator")

    print(f"decade 1e{DEC-1} -> 1e{DEC}, {len(common)} roots shared by all "
          f"{len(MODES)} datasets")
    print(f"point readings match out/summary_b15_d10.log for all "
          f"{len(MODES)} datasets\n")

    print("per-root reading, Pearson r between modes on the shared root list")
    print("  (pairing cancels variance only if these are far from zero)")
    for a, b in itertools.combinations(MODES, 2):
        xs = [data[a][r][-1] for r in common]
        ys = [data[b][r][-1] for r in common]
        print(f"    {a:<8} x {b:<8} r = {statistics.correlation(xs, ys):+.3f}")

    rng = random.Random(SEED)
    boots = {m: [] for m in MODES}
    band, sep, cons = [], [], []
    for _ in range(B):
        idx = [rng.randrange(len(common)) for _ in range(len(common))]
        p = {m: point(m, idx) for m in MODES}
        if any(v is None for v in p.values()):
            continue
        for m in MODES:
            boots[m].append(p[m])
        w = [p[m] for m in WIDE]
        band.append(max(w) - min(w))
        sep.append(p["cycq505"] - p["arith"])
        cons.append(min(p[m] for m in HIGH) - p["arith"])

    def ci(v):
        s = sorted(v)
        return s[int(0.025 * len(s))], s[int(0.975 * len(s)) - 1]

    print(f"\nmarginal bootstrap sd, then EVERY pair: the paired-difference sd")
    print("against the value it would take if the two sides were independent.")
    print(f"All {len(MODES)*(len(MODES)-1)//2}, because round 21 (R21-01) found this block printing")
    print("three, all of them pairs where nothing narrows, under a paper sentence")
    print("that generalized to every pair.")
    sd = {m: statistics.stdev(boots[m]) for m in MODES}
    for m in MODES:
        print(f"    {m:<8} sd {sd[m]:.5f}")
    print(f"\n    {'pair':<20}{'rho':>8}{'paired':>10}{'independent':>13}{'change':>9}")
    worst = 0.0
    for a, b in itertools.combinations(MODES, 2):
        rho = statistics.correlation([data[a][r][-1] for r in common],
                                     [data[b][r][-1] for r in common])
        dd = statistics.stdev([x - y for x, y in zip(boots[a], boots[b])])
        ind = math.hypot(sd[a], sd[b])
        pct = (ind - dd) / ind * 100.0
        worst = max(worst, abs(pct))
        print(f"    {a+'-'+b:<20}{rho:>+8.3f}{dd:>10.5f}{ind:>13.5f}{pct:>+8.1f}%")
    print(f"\n    joint resampling moves a pairwise sd by at most {worst:.1f}% either way.")
    sep_d = statistics.stdev([x - y for x, y in zip(boots["cycq505"], boots["arith"])])
    sep_i = math.hypot(sd["cycq505"], sd["arith"])
    print(f"    for the separation quoted in the paper (cycq505-arith) it goes")
    print(f"    {sep_i:.5f} -> {sep_d:.5f}, i.e. it WIDENS. Pairing does not help there.")

    # DIFFERENCE IN DIFFERENCES, suggested by the external referee of round 24.
    # Each structure exists at both targets, so the shift the estimator reports
    # when the target moves by 0.678 - 0.650919 can be compared against that
    # move itself, twice, in two structurally different processes.
    TARGET_SHIFT = 0.678 - 0.650919
    print(f"\ntarget shift 0.678 - 0.650919 = {TARGET_SHIFT:.6f}")
    print("observed shift minus target shift, within each matched structural pair:")
    for lowm, highm, label in (("cycq500", "cycq505", "congruence-matched"),
                               ("iidq500", "iidq505", "unconstrained")):
        d = [(x - y) - TARGET_SHIFT for x, y in zip(boots[highm], boots[lowm])]
        dpt = (pts[highm] - pts[lowm]) - TARGET_SHIFT
        dlo, dhi = ci(d)
        covers = "covers zero" if dlo <= 0 <= dhi else "EXCLUDES zero"
        print(f"    {label:<20} observed {pts[highm]-pts[lowm]:.5f}, "
              f"minus target shift {dpt:+.6f}   95% [{dlo:+.6f},{dhi:+.6f}]  {covers}")
    print("    A target move of 0.027081 therefore transmits to the estimator")
    print("    close to one for one, in both structures.")

    lo, hi = ci(sep)
    print(f"\nseparation cycq505 - arith : {pts['cycq505']-pts['arith']:.5f}"
          f"   95% [{lo:.5f}, {hi:.5f}]   (one construction, secondary)")
    clo, chi = ci(cons)
    cpt = min(pts[m] for m in HIGH) - pts["arith"]
    print(f"CONSERVATIVE separation   : {cpt:.5f}"
          f"   95% [{clo:.5f}, {chi:.5f}]")
    print(f"  = min over the two 0.678 constructions, minus the arithmetic tree,")
    print(f"    taken inside each resample rather than on the point readings.")
    print(f"    This is the quantity the paper now headlines, so it is the one")
    print(f"    that carries the interval (round 22, external critique).")
    blo, bhi = ci(band)
    w0 = max(pts[m] for m in WIDE) - min(pts[m] for m in WIDE)
    print(f"wider band width           : {w0:.5f}   95% [{blo:.5f}, {bhi:.5f}]")
    ratios = sorted((pts['cycq505'] - pts['arith']) / w for w in band if w > 1e-9)
    print(f"ratio, separation / band   : {(pts['cycq505']-pts['arith'])/w0:.1f}"
          f"   95% [{ratios[int(0.025*len(ratios))]:.1f}, "
          f"{ratios[int(0.975*len(ratios))-1]:.1f}]")
    print("\nThe separation is tight. The ratio is not a two-figure quantity, which")
    print("is why the paper quotes the separation's interval and not the ratio.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
