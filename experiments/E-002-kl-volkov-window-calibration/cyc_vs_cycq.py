#!/usr/bin/env python3
"""
E-002, part 6: is `cycq` at qval = q the same process as `cyc`?

It has to be. `cyc` runs the integer recursion w = floor((2^a u - 1)/q) with
the sibling residues advancing by c; `cycq` runs the same branching with the
value carried in log10 and divided by a real qval. At qval = q the two are
the same process up to the integer floor, which shifts a whole sibling chain
by a constant relative factor 1 - e/(2^a0 u) and is therefore negligible
away from the smallest values.

A single pair of runs showed the window estimator at 0.62829 (`cyc`) against
0.63789 (`cycq`), a gap of 0.0096. Two candidate explanations: a real
systematic inside the instrument, or Monte Carlo noise between two
independent realizations. This script settles it by running both modes at
several seeds and comparing the gap against the seed-to-seed spread of each
mode on its own.

Also reports the value-floor variant, which tests whether the real-valued
walk descending below 1 (where the arithmetic recursion cannot go) matters.

Run: python3 cyc_vs_cycq.py
"""
import math
import statistics
import subprocess

HERE = __file__.rsplit("/", 1)[0]
SEEDS = [2026, 7, 19, 44, 101, 555]


def est(args, out):
    subprocess.run([HERE + "/tree_counts", "--q", "5", "--roots", "300",
                    "--cp", "4", "8", "--buf", "9", "13", "--out", out] + args,
                   check=True, capture_output=True)
    sl = []
    for line in open(out):
        if line.startswith("#"):
            continue
        v = [int(x) for x in line.split()[2:]]
        a, b = v[1 * 5 + 4], v[4 * 5 + 4]      # N(1e5), N(1e8) at buffer 1e13
        if a > 0 and b > 0:
            sl.append(math.log10(b / a) / 3)
    return statistics.mean(sl)


def main():
    runs = {
        "cyc": ["--cyc"],
        "cycq(5.0)": ["--cycq", "5.0"],
        "cycq(5.0) floored at 1": ["--cycq", "5.0", "--minlog", "0"],
    }
    print("window estimator at buffer 1e13, 300 roots, by seed\n")
    print(f"{'mode':>24} " + " ".join(f"{s:>8}" for s in SEEDS) + f" {'mean':>8} {'sd':>7}")
    res = {}
    for name, args in runs.items():
        vals = [est(args + ["--seed", str(s)], f"{HERE}/data/tmp_cvq.txt") for s in SEEDS]
        res[name] = vals
        print(f"{name:>24} " + " ".join(f"{v:8.5f}" for v in vals)
              + f" {statistics.mean(vals):8.5f} {statistics.pstdev(vals):7.5f}")

    a, b = res["cyc"], res["cycq(5.0)"]
    gap = statistics.mean(b) - statistics.mean(a)
    sd = math.sqrt(statistics.pvariance(a) + statistics.pvariance(b)) / math.sqrt(len(SEEDS))
    print(f"\ncycq(5.0) minus cyc: {gap:+.5f}, standard error of that difference {sd:.5f}"
          f"  ({abs(gap)/sd:.1f} sigma)")
    print("The single-pair gap that prompted this check was +0.0096.")


if __name__ == "__main__":
    main()
