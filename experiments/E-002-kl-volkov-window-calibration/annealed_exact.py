#!/usr/bin/env python3
"""
E-002, part 1: the exact annealed counting function of the i.i.d. model,
in closed form, with no simulation.

For the branching random walk whose annealed pressure is
rho(alpha) = q^(alpha-1)/(2^alpha-1), a level-k node reached by exponents
a_1..a_k (each a_i >= 1) sits at value  u0 * 2^A / q^k  with A = sum a_i,
and the expected number of children of a fertile node at any prescribed
exponent n >= 1 is exactly 1/q. Hence

    E[# level-k nodes with sum a_i = A] = q^(-k) * C(A-1, k-1).

Summing over A <= N_k(t) := floor((t + k log10 q)/log10 2) and applying the
hockey-stick identity sum_{A=k}^{N} C(A-1,k-1) = C(N,k) gives

    M(t) := E[ N(u0 * 10^t) ] = sum_{k>=1} C(N_k(t), k) / q^k .

Exact, one line, and evaluable at t far beyond any enumeration. This is the
instrument used to decide whether the Kontorovich-Lagarias versus Volkov
dispute is reachable by direct counting at all.

Run: python3 annealed_exact.py [q]
"""
import math
import sys

L2 = math.log10(2.0)


def alpha_minus(q):
    """smaller root of q^(alpha-1) = 2^alpha - 1"""
    def f(a):
        return (q ** (a - 1.0)) - (2.0 ** a - 1.0)
    lo, hi = 1e-9, 1.0 - 1e-12
    # f(lo) > 0, f(1) = 0; bisect on the crossing strictly below 1
    if q == 3:
        return 1.0
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def logM(t, q, kmax=None, qval=None):
    """log10 of the exact annealed count at value ratio 10^t. Log-space sum.

    `q` sets the offspring intensity 1/q per exponent; `qval` sets the value
    denominator and defaults to q. Mode `cycq` of the enumerator separates
    them, which is how its counting exponent is tuned: the exponent solves
    qval^alpha = q (2^alpha - 1).
    """
    # log10 of the VALUE denominator, which sets how deep a level reaches;
    # the per-level weight is q^(-k) and uses the intensity denominator, a
    # distinction that only shows up when the two differ
    lq = math.log10(q if qval is None else qval)
    lqw = math.log10(q)
    if kmax is None:
        # mass per level decays like ~0.978^k at q=5, so go deep; checked
        # against kmax = 8000 at t = 1, 3, 8, 20, agreeing to 8 decimals
        kmax = max(1500, int(400 + 60 * t))
    terms = []
    for k in range(1, kmax + 1):
        N = int(math.floor((t + k * lq) / L2 + 1e-12))
        if N < k:
            continue
        # log10 C(N,k) - k log10 q
        lc = (math.lgamma(N + 1) - math.lgamma(k + 1) - math.lgamma(N - k + 1)) / math.log(10.0)
        terms.append(lc - k * lqw)
    if not terms:
        return None
    mx = max(terms)
    s = sum(10.0 ** (x - mx) for x in terms)
    return mx + math.log10(s)


def main():
    q = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    am = alpha_minus(q)
    print(f"q={q}   alpha_minus = {am:.6f}   (KL eta_{q},BP)")
    print(f"exact annealed local slope  d log10 M / d log10 t, per decade of x/u0\n")
    print(f"{'t=log10(x/u0)':>14} {'slope on [t,t+1]':>18} {'deficit vs alpha-':>18}")

    ts = [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20, 30, 50, 80, 120, 200, 350, 600, 1000]
    prev = {}
    for t in ts:
        a = logM(t, q)
        b = logM(t + 1, q)
        if a is None or b is None:
            continue
        sl = b - a
        prev[t] = sl
        print(f"{t:>14} {sl:>18.6f} {am - sl:>18.6f}")

    # how deep must the window be for the local slope to sit within
    # half the KL/Volkov separation of alpha_minus?
    if q == 5:
        gap = 0.678 - 0.650919
        print(f"\nKL vs Volkov separation Delta = {gap:.6f}")
        for tol_name, tol in (("Delta", gap), ("Delta/2", gap / 2), ("Delta/4", gap / 4)):
            t = 1.0
            found = None
            while t < 200000:
                sl = logM(t + 1, q) - logM(t, q)
                if am - sl < tol:
                    found = t
                    break
                t *= 1.35
            if found:
                print(f"  local slope comes within {tol_name} = {tol:.5f} of alpha_- "
                      f"at t ~ {found:.0f}, i.e. x/u0 ~ 10^{found:.0f}")
            else:
                print(f"  local slope never comes within {tol_name} below t = 200000")


if __name__ == "__main__":
    main()
