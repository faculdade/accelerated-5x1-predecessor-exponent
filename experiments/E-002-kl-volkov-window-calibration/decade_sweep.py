#!/usr/bin/env python3
"""
E-002, part 7: the separation at every checkpoint decade the b15 grid holds.

The paper cross-checked its calibration decade against one shallower decade and
concluded that the separation does not shrink with depth. The grid holds five,
and round 19 (R19-03) found that the two unreported readings at the cross-check
decade are the two that set the band, so the unit quoted there is 4.3x wider
than at the calibration decade, and the gap sequence is not monotone.

Round 20 then found that the two claims round 19 wrote from this table are both
false, so the script now computes what it takes to state them properly.

  R20-02: "the separation is present at every decade" fails at the shallowest.
  There the arithmetic reading sits ABOVE the whole calibration band rather than
  inside it, an inversion a table of gaps alone hides, and the gap's own 95%
  interval covers zero.

  R20-03: "the band narrows monotonically" holds with probability about 0.48
  under this same bootstrap, a coin flip. The two deepest steps are solid
  (0.996, 0.990); the first is not (0.62).

So the table now prints, per decade, whether the arithmetic reading is inside
the band and the gap's bootstrap interval, and it prints the monotonicity
probability step by step instead of asserting the ordering.

Estimator: `summary.py`'s own `decade_estimate`, imported. The two decades the
paper already quotes must match `out/summary_b15_d8.log` and
`out/summary_b15_d10.log`; the script checks that and aborts if not.

Serial (Rule 9c): five decades x seven datasets, plus 2000 joint resamples per
decade, over committed count files; about 7 s, under the one-minute line.

Run: python3 decade_sweep.py
"""
import random
import sys

sys.path.insert(0, ".")
from summary import decade_estimate, load, slope, aitken, mean   # noqa: E402

MODES = ["cycq500", "cyc", "iid", "iidq500", "arith", "cycq505", "iidq505"]
WIDE = ["cycq500", "cyc", "iid", "iidq500"]      # all four built to 0.650919
HIGH = ["cycq505", "iidq505"]                    # both built to 0.678      # the three built to 0.650919
CHECK = {(8, "arith"): 0.64622, (8, "cycq505"): 0.67200,
         (10, "arith"): 0.64926, (10, "cycq505"): 0.67748}


def curves(mode, dec):
    _, roots, mats, cp_lo, buf_lo, n_cp, n_buf, _ = load(f"data/q5_{mode}_b15.txt")
    ci = dec - 1 - cp_lo
    use = [bi for bi in range(n_buf) if buf_lo + bi > dec]
    out = {}
    for r, m in zip(roots, mats):
        v = [slope(m, ci, ci + 1, bi) for bi in use]
        if all(x is not None for x in v):
            out[r] = v
    return out, len(use)


def main():
    print("the arithmetic tree against the construction retuned to 0.678, at")
    print("every checkpoint decade the b15 grid supports\n")
    print(f"{'decade':<12}" + "".join(f"{m:>10}" for m in MODES)
          + f"{'band':>9}{'to floor':>10}  arith in band?  to-floor 95%")
    bands_boot, bands_pt = [], []
    for dec in (6, 7, 8, 9, 10):
        pts = {}
        for m in MODES:
            p, _, _ = decade_estimate(f"data/q5_{m}_b15.txt", dec_top=dec)
            pts[m] = p
            if (dec, m) in CHECK and round(p, 5) != CHECK[(dec, m)]:
                raise SystemExit(f"ERROR: decade {dec} {m} reads {p:.5f}, the "
                                 f"committed log says {CHECK[(dec, m)]:.5f}")
        w = [pts[m] for m in WIDE]
        band, top, bot = max(w) - min(w), max(w), min(w)
        # conservative: to the FLOOR of the 0.678 family, the smallest gap the
        # two admit. Until round 24 this line used cycq505 alone, so this
        # script and referee_tables.py answered the same question differently.
        gap = min(pts[m] for m in HIGH) - pts["arith"]
        inside = bot <= pts["arith"] <= top

        d, nb = {}, None
        for m in MODES:
            d[m], nb = curves(m, dec)
        common = sorted(set.intersection(*(set(v) for v in d.values())))

        def pt(m, idx):
            rows = [d[m][common[i]] for i in idx]
            return aitken([mean([r[b] for r in rows]) for b in range(nb)])

        rng = random.Random(17)
        gaps, bb = [], []
        for _ in range(2000):
            idx = [rng.randrange(len(common)) for _ in common]
            q = {m: pt(m, idx) for m in MODES}
            if any(v is None for v in q.values()):
                continue
            gaps.append(min(q[m] for m in HIGH) - q["arith"])
            bb.append(max(q[m] for m in WIDE) - min(q[m] for m in WIDE))
        gaps.sort()
        lo, hi = gaps[int(0.025 * len(gaps))], gaps[int(0.975 * len(gaps)) - 1]
        bands_boot.append(bb)
        bands_pt.append(band)
        flag = "yes" if inside else "NO, ABOVE" if pts["arith"] > top else "NO, BELOW"
        print(f"1e{dec-1}->1e{dec:<7}" + "".join(f"{pts[m]:>10.5f}" for m in MODES)
              + f"{band:>9.5f}{gap:>9.4f}  {flag:<14} [{lo:+.5f},{hi:+.5f}]")

    print("\nAt the shallowest decade the arithmetic reading sits ABOVE the")
    print("0.650919 span and INSIDE the 0.678 span, so the conservative gap is")
    print("negative there and its interval is entirely negative; that decade")
    print("carries no calibration weight. At 1e6->1e7 the interval still covers")
    print("zero. It excludes zero in the positive direction only from 1e7->1e8")
    print("onward. Until round 25 this block asserted the opposite of both, a")
    print("narrative written for the pre-conservative gap and never re-read")
    print("against the numbers printed directly above it.")

    n = min(len(b) for b in bands_boot)
    allmono = sum(1 for i in range(n)
                  if all(bands_boot[j][i] > bands_boot[j + 1][i] for j in range(4)))
    print(f"\nP(the band narrows at every step) = {allmono/n:.3f}, so the ordering")
    print("is not a fact to assert. Step by step:")
    for j in range(4):
        s = sum(1 for i in range(n) if bands_boot[j][i] > bands_boot[j + 1][i])
        print(f"   step {j+1}  {bands_pt[j]:.5f} -> {bands_pt[j+1]:.5f}   P = {s/n:.3f}")
    print(f"\nWhat IS solid: the band at 1e6->1e7 is {bands_pt[1]:.5f} and at")
    print(f"1e9->1e10 is {bands_pt[4]:.5f}, a factor of {bands_pt[1]/bands_pt[4]:.1f}.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
