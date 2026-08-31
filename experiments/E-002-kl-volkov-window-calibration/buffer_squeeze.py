#!/usr/bin/env python3
"""
E-002, part 5: how much does the truncation extrapolation itself cost?

The deepest decades of the deep run have the fewest buffers above them, so
their Aitken step is the least constrained: at decade 1e11->1e12 with
buffers to 1e17, the three points feeding Aitken sit only 3, 4 and 5
decades above the decade's top. The bootstrap band on those decades
measures root resampling and says nothing about this.

So take a decade that has nine buffers, redo its Aitken using only the
three buffers that reproduce the deep decades' headroom, and compare. The
difference is the extrapolation error the deep decades are exposed to.

Run: python3 buffer_squeeze.py data/q5_arith_b17.txt
"""
import sys

from analyze import aitken, load, mean, slope


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "data/q5_arith_b17.txt"
    hdr, roots, mats, cp_lo, buf_lo, n_cp, n_buf, _ = load(path)
    print(f"=== buffer-squeeze test on {path} (mode {hdr['mode']}) ===")
    print("For each decade: Aitken on every usable buffer, versus Aitken on only")
    print("the three buffers sitting 3, 4 and 5 decades above the decade's top,")
    print("which is all the deepest decade of this run has.\n")
    print(f"{'decade':>14} {'all buffers':>12} {'headroom 3-5':>13} {'difference':>11}")
    for ci in range(n_cp - 1):
        top = cp_lo + ci + 1
        row = []
        for bi in range(n_buf):
            v = [slope(m, ci, ci + 1, bi) for m in mats]
            row.append(mean(v) if all(x is not None for x in v) else None)
        full = [row[bi] for bi in range(n_buf)
                if row[bi] is not None and buf_lo + bi > top]
        squeeze = [row[bi] for bi in range(n_buf)
                   if row[bi] is not None and top + 3 <= buf_lo + bi <= top + 5]
        if len(full) < 3 or len(squeeze) < 3:
            continue
        a, b = aitken(full), aitken(squeeze)
        print(f"  1e{top-1}->1e{top:<8} {a:12.4f} {b:13.4f} {b-a:+11.4f}")
    print("\nA positive difference means the squeezed extrapolation over-reads, so")
    print("the deep decades of this run carry that much upward extrapolation error")
    print("on top of their bootstrap band.")


if __name__ == "__main__":
    main()
