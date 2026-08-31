#!/usr/bin/env python3
"""
E-002, part 2: check that the i.i.d. simulator reproduces the exact annealed
count of `annealed_exact.py`.

One correction has to be applied before they can be compared. The closed
form lets the root be sterile with probability 1/q, while the enumerator
draws a fertile root residue in every mode, to match the arithmetic run
whose roots are sampled with u mod q != 0. Conditioning the root to be
fertile raises the level-1 intensity from 1/q to 1/d = 1/(q-1) per
exponent, and every deeper level inherits the same factor, so the
simulator's mean is q/(q-1) times the closed form. That factor is 1.25 at
q = 5 and is divided out below.

After that, WHAT THE COMPARISON MEANS DEPENDS ON WHICH MODE IS COMPARED, and
until round 27 this docstring collapsed the two cases into one claim that they
"must agree".

For the `cycq` modes the closed form is EXACT: those are the multiplicative
relaxations the derivation is written for, with the value denominator qval, so
agreement up to sampling error is a genuine validation of the closed form.

For the `--iid` integer mode it is NOT an exact identity. That process seeds
each sibling family with floor((2^a0 u - 1)/q) and continues it by the affine
recurrence w_{j+1} = 2^d w_j + (2^d-1)/q, which is a perturbation of the
multiplicative walk V_k = u0 2^A / q^k the closed form describes; its children
are never larger than the multiplicative ones. So for `--iid` this is a
comparison against the multiplicative annealed baseline, not an equality that
must hold, and the ratios below sitting under 1 are consistent with both the
perturbation and the sampling error described next.

The sampling error is real either way:
the total progeny of this branching random walk has tail index
alpha_+/alpha_- = 1/0.650919 = 1.5363, so its variance is infinite and the
finite-sample mean is highly unstable: it can sit substantially above or
below the true mean and moves by several percent from seed to seed. The
table below shows both directions, seed 33 landing above the closed form at
every t and seed 22 below. An earlier version of this note claimed the
sample mean always sits BELOW, which two of its own four seeds refute. The comparison also shows the truncation dependence one would expect: raising
the buffer moves the simulated mean toward the multiplicative baseline, from
0.7240 at 1e11 to 0.9191 at 1e13 for seed 44 at t=4. An earlier version of
this note claimed the deficit "does not depend on the truncation buffer",
which those three lines refute directly. What the check is for is narrower:
the absence of a discrepancy growing with the sample size is a useful sanity
check, but not a proof of implementation correctness. For the integer iid
mode this stays a sanity comparison against the multiplicative baseline, not
an exact identity.

Run: python3 check_mean_vs_annealed.py
"""
import math
import statistics
import subprocess
import sys

import annealed_exact as ae

HERE = __file__.rsplit("/", 1)[0]
ROOT = 1000003
FERTILE_ROOT = 5.0 / 4.0   # q/(q-1): the enumerator always starts fertile


def run(seed, nroots, cp_lo, cp_hi, buf_lo, buf_hi, out):
    subprocess.run([HERE + "/tree_counts", "--q", "5", "--iid",
                    "--fixedroot", str(ROOT), "--roots", str(nroots),
                    "--seed", str(seed), "--cp", str(cp_lo), str(cp_hi),
                    "--buf", str(buf_lo), str(buf_hi), "--out", out],
                   check=True, capture_output=True)
    rows = []
    for line in open(out):
        if line.startswith("#"):
            continue
        rows.append([int(x) for x in line.split()[2:]])
    return rows


def main():
    n_buf = 3
    print("i.i.d. simulator mean count vs exact annealed M(t), root =", ROOT)
    print(f"{'seed':>6} {'n':>8} " + "  ".join(f"t={t}" for t in (1, 2, 3, 4)))
    for seed, n in ((11, 20000), (22, 20000), (33, 20000), (44, 200000)):
        rows = run(seed, n, 7, 10, 11, 13, f"{HERE}/data/mc_{seed}_{n}.txt")
        cells = []
        for ci in range(4):
            t = math.log10(10 ** (7 + ci) / ROOT)
            m = statistics.mean(r[ci * n_buf + 2] for r in rows) - 1.0
            ex = FERTILE_ROOT * 10 ** ae.logM(t, 5)
            cells.append(f"{m/ex:6.4f}")
        print(f"{seed:>6} {n:>8} " + "  ".join(cells))

    print("\nsame quantity at three truncation buffers (seed 44, n=200000):")
    rows = run(44, 200000, 7, 10, 11, 13, f"{HERE}/data/mc_44_200000.txt")
    for bi, b in enumerate((11, 12, 13)):
        cells = []
        for ci in range(4):
            t = math.log10(10 ** (7 + ci) / ROOT)
            m = statistics.mean(r[ci * n_buf + bi] for r in rows) - 1.0
            ex = FERTILE_ROOT * 10 ** ae.logM(t, 5)
            cells.append(f"{m/ex:6.4f}")
        print(f"  buffer 1e{b}: " + "  ".join(cells))

    # the tunable-exponent mode, against the same closed form with qval
    print("\nmode cycq against the closed form with a separate value denominator:")
    for qv in ("5.00000", "5.05398"):
        subprocess.run([HERE + "/tree_counts", "--q", "5", "--cycq", qv,
                        "--fixedroot", str(ROOT), "--roots", "40000",
                        "--seed", "77", "--cp", "7", "10", "--buf", "11", "13",
                        "--out", f"{HERE}/data/mc_cycq_{qv}.txt"],
                       check=True, capture_output=True)
        rows = []
        for line in open(f"{HERE}/data/mc_cycq_{qv}.txt"):
            if not line.startswith("#"):
                rows.append([int(x) for x in line.split()[2:]])
        cells = []
        for ci in range(4):
            t = math.log10(10 ** (7 + ci) / ROOT)
            m = statistics.mean(r[ci * n_buf + 2] for r in rows) - 1.0
            ex = FERTILE_ROOT * 10 ** ae.logM(t, 5, qval=float(qv))
            cells.append(f"{m/ex:6.4f}")
        print(f"  qval={qv}: " + "  ".join(cells))


if __name__ == "__main__":
    sys.exit(main())
