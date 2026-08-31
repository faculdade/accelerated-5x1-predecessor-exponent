#!/usr/bin/env python3
"""
E-002, part 4: is the arithmetic tree's low spread a real property, or just
bookkeeping?

Mode `arith` gives one deterministic tree per root, while `cyc` and `iid`
give one draw from an ensemble per root. So the controls carry realization
noise that the arithmetic tree cannot carry, and a smaller spread in `arith`
proves nothing by itself.

This compares like with like. For a handful of the actual arithmetic roots,
it draws many control realizations from that same root and measures

  - the within-root spread of log10 N(1e8) under each control;
  - where the one arithmetic value falls inside that control distribution.

If the within-root spread alone is as large as the across-root spread of the
arithmetic run, the variance-reduction reading dissolves. If the arithmetic
value also sits systematically high in the control distribution, the
arithmetic tree is not a typical realization of its own branching model.

Run: python3 within_root_spread.py
"""
import math
import statistics
import subprocess

HERE = __file__.rsplit("/", 1)[0]
NREAL = 3000
# a shallower truncation buffer than the main runs: this compares spreads,
# and every mode here is truncated identically, so the comparison holds
CP_LO, CP_HI, BUF_LO, BUF_HI = 4, 8, 9, 11
N_BUF = BUF_HI - BUF_LO + 1
I5, I8 = 5 - CP_LO, 8 - CP_LO


def counts(args, out):
    subprocess.run([HERE + "/tree_counts", "--q", "5", "--cp", str(CP_LO), str(CP_HI),
                    "--buf", str(BUF_LO), str(BUF_HI), "--out", out] + args,
                   check=True, capture_output=True)
    rows = []
    for line in open(out):
        if line.startswith("#"):
            continue
        v = [int(x) for x in line.split()[2:]]
        rows.append((v[I5 * N_BUF + N_BUF - 1], v[I8 * N_BUF + N_BUF - 1]))
    return rows


def main():
    # the first roots of the standard arith sample
    arith_all = counts(["--roots", "300"], f"{HERE}/data/tmp_arith.txt")
    roots = []
    for line in open(f"{HERE}/data/tmp_arith.txt"):
        if not line.startswith("#"):
            roots.append(int(line.split()[0]))
    sd_across = statistics.pstdev([math.log10(b) for _, b in arith_all if b > 0])
    sl_across = statistics.pstdev([math.log10(b / a) / 3 for a, b in arith_all if a > 0 and b > 0])
    print(f"arith, across 300 deterministic roots: "
          f"sd log10 N(1e8) = {sd_across:.4f}, sd of per-root slope = {sl_across:.4f}\n")

    print(f"{'root':>7} {'arith log10N':>13} {'arith slope':>12}   "
          f"{'ctrl':>5} {'within sd log10N':>17} {'within sd slope':>16} "
          f"{'pct of arith N':>15} {'pct of arith slope':>19}")
    sds_n, sds_s, pn, ps = {}, {}, {}, {}
    for k, root in enumerate(roots[:10]):
        a5, a8 = arith_all[k]
        av = math.log10(a8) if a8 > 0 else None
        asl = math.log10(a8 / a5) / 3 if a5 > 0 and a8 > 0 else None
        for mode, flag in (("cyc", "--cyc"), ("iid", "--iid")):
            rows = counts([flag, "--fixedroot", str(root), "--roots", str(NREAL),
                           "--seed", str(9000 + k)], f"{HERE}/data/tmp_{mode}_{root}.txt")
            lg = [math.log10(b) for _, b in rows if b > 0]
            sl = [math.log10(b / a) / 3 for a, b in rows if a > 0 and b > 0]
            sd_n, sd_s = statistics.pstdev(lg), statistics.pstdev(sl)
            p_n = 100.0 * sum(1 for v in lg if v < av) / len(lg)
            p_s = 100.0 * sum(1 for v in sl if v < asl) / len(sl)
            sds_n.setdefault(mode, []).append(sd_n)
            sds_s.setdefault(mode, []).append(sd_s)
            pn.setdefault(mode, []).append(p_n)
            ps.setdefault(mode, []).append(p_s)
            head = f"{root:>7} {av:>13.4f} {asl:>12.5f}" if mode == "cyc" else " " * 34
            print(f"{head}   {mode:>5} {sd_n:>17.4f} {sd_s:>16.5f} "
                  f"{p_n:>14.1f}% {p_s:>18.1f}%")

    print()
    for mode in ("cyc", "iid"):
        print(f"{mode}: mean within-root sd log10 N(1e8) = {statistics.mean(sds_n[mode]):.4f}, "
              f"mean within-root sd of slope = {statistics.mean(sds_s[mode]):.4f}, "
              f"median percentile of the arithmetic value = "
              f"{statistics.median(pn[mode]):.1f}% (count), {statistics.median(ps[mode]):.1f}% (slope)")


if __name__ == "__main__":
    main()
