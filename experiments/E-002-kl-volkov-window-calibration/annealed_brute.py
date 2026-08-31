#!/usr/bin/env python3
"""
E-002, part 1b: the brute-force double sum behind `annealed_exact.py`.

`annealed_exact.py` evaluates

    M(t) = sum_{k>=1} C(N_k(t), k) / q^k ,   N_k(t) = floor((t + k log10 q)/log10 2),

which it reaches from the double sum

    M(t) = sum_{k>=1} q^(-k) sum_{A=k}^{N_k(t)} #{(a_1..a_k) : a_i >= 1, sum a_i = A}

by three steps: the admission bound 2^A/q^k <= 10^t becomes A <= N_k(t); the
count of exponent tuples is C(A-1, k-1); and the hockey-stick identity collapses
sum_{A=k}^{N} C(A-1,k-1) into C(N,k). This script checks ALL THREE against
routes that use none of them.

Until round 17 it checked only the third. Both sides took the SAME N from the
same floating-point `n_k`, so a +/-1 error in the bound cancelled and the check
still reported MATCH (R17-03). The brute side now derives its own bound in exact
integer arithmetic, `max{A : 2^A <= 10^t q^k}`, straight from the admission
inequality, with no logarithm anywhere; `n_k`'s float result is compared against
it rather than trusted.

The inner count is rebuilt by dynamic programming over the tuples themselves,

    f[1][a] = 1  (a >= 1),      f[k][A] = sum_{a>=1} f[k-1][A-a],

which is the definition, not the formula. Everything runs in `fractions.Fraction`,
so the comparison is exact: an identity checked in floating point is a claim about
rounding, not about the identity (CLAUDE.md Rule 11c).

Two checks per (q, t): each level k separately, and the truncated total.

Run: python3 annealed_brute.py            # q = 3, 5, 7 and t = 1..4
     python3 annealed_brute.py --selftest # negative control, must report FAIL
"""
import math
import sys
from fractions import Fraction

L2 = math.log10(2.0)
KMAX = 60


def n_k(t, k, q):
    """floor((t + k log10 q)/log10 2), the deepest exponent sum a level-k node reaches.

    Floating point, with a fudge term; this is the form `annealed_exact.py` uses
    and therefore the form under test."""
    return int(math.floor((t + k * math.log10(q)) / L2 + 1e-12))


def n_k_exact(t, k, q):
    """the same bound from the admission inequality itself, in exact integers.

    A level-k node with exponent sum A is admitted when 2^A / q^k <= 10^t, i.e.
    2^A <= 10^t q^k. For integer t and k that right-hand side is an exact
    integer, so the largest admissible A is its bit_length minus one. No
    logarithm, no rounding, nothing shared with `n_k`."""
    return (10 ** t * q ** k).bit_length() - 1


def tuple_counts(kmax, amax):
    """f[k][A] = #{(a_1..a_k): a_i >= 1, sum a_i = A}, by DP over the tuples.

    Uses no binomial coefficient and no hockey-stick identity; this is the
    independent route the closed form is checked against.
    """
    f = [[0] * (amax + 1) for _ in range(kmax + 1)]
    for a in range(1, amax + 1):
        f[1][a] = 1
    for k in range(2, kmax + 1):
        for A in range(k, amax + 1):
            f[k][A] = sum(f[k - 1][A - a] for a in range(1, A - k + 2))
    return f


def report(q, t, f, bug=None):
    """Compare both routes at (q, t).

    `bug` is the negative control and names WHICH link to break:
      "bound"  perturbs the closed form's admission bound N_k;
      "count"  perturbs the tuple count the DP supplies to the brute side;
      "hockey" perturbs C(N,k).
    Each link needs its own probe: a control that only breaks one link cannot
    show that the others are checked at all (R17-03). Round 18 (R18-07) found
    this header announcing one probe per link while shipping two for three."""
    amax = n_k_exact(t, KMAX, q)
    brute = Fraction(0)
    closed = Fraction(0)
    bad_levels = []
    bad_bounds = []
    for k in range(1, KMAX + 1):
        N_ex = n_k_exact(t, k, q)                 # brute side: from the inequality
        N_fl = n_k(t, k, q)                       # closed side: the form under test
        if bug == "bound" and k == 3:
            N_fl += 1
        if N_fl != N_ex:
            bad_bounds.append((k, N_ex, N_fl))
        if N_ex < k and N_fl < k:
            continue
        inner = sum(f[k][A] for A in range(k, N_ex + 1))       # double sum, own bound
        if bug == "count" and k == 3:
            inner += 1
        hock = math.comb(N_fl + (1 if bug == "hockey" and k == 3 else 0), k)
        if inner != hock:
            bad_levels.append((k, inner, hock))
        brute += Fraction(inner, q ** k)
        closed += Fraction(hock, q ** k)
    ok = (not bad_levels) and (not bad_bounds) and (brute == closed)
    why = ""
    if bad_bounds:
        why = (f"  bound differs at level {bad_bounds[0][0]}: exact "
               f"{bad_bounds[0][1]} vs n_k {bad_bounds[0][2]}")
    elif bad_levels:
        why = (f"  first bad level {bad_levels[0][0]}: "
               f"{bad_levels[0][1]} vs {bad_levels[0][2]}")
    print(f"  q={q} t={t}: levels 1..{KMAX}, A up to {amax}, "
          f"{'MATCH' if ok else 'FAIL'} exactly" + why)
    if ok:
        print(f"            truncated M(t) = {float(brute):.12f} by both routes")
    return ok


def main():
    qs, ts = [3, 5, 7], [1, 2, 3, 4]
    amax = max(n_k_exact(max(ts), KMAX, q) for q in qs)
    f = tuple_counts(KMAX, amax)

    if "--selftest" in sys.argv:
        # one probe per link. A control that breaks only the hockey-stick step
        # says nothing about whether the bound is checked, which is exactly how
        # the first version of this script passed while blind to it (R17-03).
        print("NEGATIVE CONTROLS, one per link of the chain\n")
        fail = 0
        for link, what in (("bound", "the admission bound N_k, perturbed at k=3"),
                           ("count", "the DP's tuple count, perturbed at k=3"),
                           ("hockey", "the hockey-stick collapse C(N,k), perturbed at k=3")):
            allok = True
            for q in qs:
                for t in ts:
                    allok &= report(q, t, f, bug=link)
            verdict = "PASS, the check notices" if not allok else "BROKEN, the check is blind"
            print(f"  -> {what}: {verdict}\n")
            fail |= 1 if allok else 0
        sys.exit(fail)

    print("brute-force double sum vs closed form, exact rational arithmetic")
    print(f"DP table: k up to {KMAX}, A up to {amax}")
    print("the brute side takes its bound from 2^A <= 10^t q^k in exact integers;")
    print("the closed side takes it from the floating-point n_k under test\n")
    allok = True
    for q in qs:
        for t in ts:
            allok &= report(q, t, f)
        print()
    print("all (q,t) agree exactly" if allok else "MISMATCH")
    sys.exit(0 if allok else 1)


if __name__ == "__main__":
    main()
