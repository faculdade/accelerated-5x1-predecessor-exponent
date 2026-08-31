#!/usr/bin/env python3
"""
E-002, part 3: where the arithmetic tree and its stochastic controls differ.

All three runs use the same 300 roots, the same window and the same
truncation buffers, so the columns are directly comparable.

Run: python3 compare_modes.py data/a.txt data/b.txt ...
"""
import math
import statistics
import sys

from analyze import load


def main():
    print(f"{'mode':>7} {'mean N(1e5)':>12} {'mean N(1e8)':>12} "
          f"{'gmean N(1e5)':>13} {'gmean N(1e8)':>13} {'N=0 at 1e5':>11} "
          f"{'sd log10 N(1e8)':>16} {'mean slope':>11}")
    for path in sys.argv[1:]:
        hdr, roots, mats, cp_lo, buf_lo, n_cp, n_buf, _ = load(path)
        bi = n_buf - 1
        a, b = 5 - cp_lo, 8 - cp_lo
        n5 = [m[a][bi] for m in mats]
        n8 = [m[b][bi] for m in mats]
        zero5 = sum(1 for v in n5 if v == 0)
        nz = [(x, y) for x, y in zip(n5, n8) if x > 0 and y > 0]
        g5 = statistics.mean(math.log10(x) for x, _ in nz)
        g8 = statistics.mean(math.log10(y) for _, y in nz)
        sd8 = statistics.pstdev([math.log10(y) for _, y in nz])
        sl = statistics.mean(math.log10(y / x) / 3 for x, y in nz)
        print(f"{hdr['mode']:>7} {statistics.mean(n5):12.1f} {statistics.mean(n8):12.1f} "
              f"{g5:13.4f} {g8:13.4f} {zero5:11d} {sd8:16.4f} {sl:11.5f}")

    print("\n(gmean columns are means of log10 N, i.e. logs of geometric means)")


if __name__ == "__main__":
    main()
