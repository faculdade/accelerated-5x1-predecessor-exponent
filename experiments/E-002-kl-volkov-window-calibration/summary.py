#!/usr/bin/env python3
"""
E-002, part 7: the comparison table.

Every row is the same estimator on the same 300 roots with the same window
and the same truncation buffers. Only the process differs. Two of the rows
are processes tuned to the indicated annealed pressure target. Equality of
realized and pressure exponents is NOT assumed by this script and is not
proved by the experiment; the paper carries it as a stated assumption. Their
readings say what the estimator returns for each of the two disputed
values, and the arithmetic row can be read against them.

Run: python3 summary.py            # the b13 grid, window 1e5..1e8
     python3 summary.py b15 10     # the b15 grid, deepest decade 1e9->1e10
"""
import math
import random
import sys

from analyze import aitken, load, mean, slope

HERE = __file__.rsplit("/", 1)[0]

SUFFIX = sys.argv[1] if len(sys.argv) > 1 else "b13"
DEC_TOP = int(sys.argv[2]) if len(sys.argv) > 2 else 8

ROWS = [
    ("cycq 5.00000", f"q5_cycq500_{SUFFIX}.txt", 0.650919),
    ("cycq 5.05398", f"q5_cycq505_{SUFFIX}.txt", 0.678000),
    ("cyc",          f"q5_cyc_{SUFFIX}.txt",     0.650919),
    ("iid",          f"q5_iid_{SUFFIX}.txt",     0.650919),
    ("iidq 5.00000", f"q5_iidq500_{SUFFIX}.txt", 0.650919),
    ("iidq 5.05398", f"q5_iidq505_{SUFFIX}.txt", 0.678000),
    ("arith",        f"q5_arith_{SUFFIX}.txt",   None),
]


def window_estimate(path, seed=3):
    hdr, roots, mats, cp_lo, buf_lo, n_cp, n_buf, _ = load(path)
    a, b = 5 - cp_lo, 8 - cp_lo
    per_root = []
    for m in mats:
        r = [slope(m, a, b, bi) for bi in range(n_buf)]
        if all(v is not None for v in r):
            per_root.append(r)
    curve = [mean([r[bi] for r in per_root]) for bi in range(n_buf)]
    point = aitken(curve)
    rng = random.Random(seed)
    boots = []
    for _ in range(2000):
        idx = [rng.randrange(len(per_root)) for _ in range(len(per_root))]
        c = [mean([per_root[i][bi] for i in idx]) for bi in range(n_buf)]
        v = aitken(c)
        if v is not None:
            boots.append(v)
    boots.sort()
    return point, boots[int(0.025 * len(boots))], boots[int(0.975 * len(boots)) - 1]


def decade_estimate(path, dec_top=8, seed=17):
    hdr, roots, mats, cp_lo, buf_lo, n_cp, n_buf, _ = load(path)
    ci = dec_top - 1 - cp_lo
    per_root = []
    for m in mats:
        r = [slope(m, ci, ci + 1, bi) for bi in range(n_buf)]
        if all(v is not None for v in r):
            per_root.append(r)
    use = [bi for bi in range(n_buf) if buf_lo + bi > dec_top]
    curve = [mean([r[bi] for r in per_root]) for bi in use]
    point = aitken(curve)
    rng = random.Random(seed)
    boots = []
    for _ in range(2000):
        idx = [rng.randrange(len(per_root)) for _ in range(len(per_root))]
        c = [mean([per_root[i][bi] for i in idx]) for bi in use]
        v = aitken(c)
        if v is not None:
            boots.append(v)
    boots.sort()
    return point, boots[int(0.025 * len(boots))], boots[int(0.975 * len(boots)) - 1]


def main():
    print(f"Same estimator, same 300 roots, window 1e5..1e8, grid {SUFFIX},")
    print("Aitken in the truncation buffer, 2000-resample bootstrap over roots.\n")
    dl = f"decade 1e{DEC_TOP-1}->1e{DEC_TOP}"
    print(f"{'process':>14} {'pressure target':>16}   {'window estimator':>26}   {dl:>26}")
    for name, fn, truth in ROWS:
        try:
            w = window_estimate(f"{HERE}/data/{fn}")
            d = decade_estimate(f"{HERE}/data/{fn}", dec_top=DEC_TOP)
        except FileNotFoundError:
            print(f"{name:>14}  (missing {fn})")
            continue
        tt = f"{truth:.6f}" if truth else "disputed"
        print(f"{name:>14} {tt:>14}   {w[0]:.5f} [{w[1]:.5f},{w[2]:.5f}]   "
              f"{d[0]:.5f} [{d[1]:.5f},{d[2]:.5f}]")
    print("\nThe two cycq rows bracket what the estimator returns for the two")
    print("disputed exponents. The arithmetic row is read against them.")


if __name__ == "__main__":
    sys.exit(main())
