#!/usr/bin/env python3
"""
E-003: does the MULTITYPE (matrix) growth rate of the accelerated qx+1
reverse tree equal the scalar annealed exponent alpha_- ?

Why this exists. E-002's theory is scalar: it assumes a fertile node has
expected 1/q children at each exponent n >= 1, which gives the annealed
pressure rho(alpha) = q^(alpha-1)/(2^alpha - 1) and the Kontorovich-Lagarias
value alpha_- as its smaller root. That assumption averages over residues.
H-002 proved the residues are NOT averaged: siblings walk a deterministic
arithmetic progression mod q, so which siblings are fertile is determined,
not sampled.

Menshikov, Petritis and Volkov (Bernoulli 13(4), 2007) answer the
corresponding counting question for a coloured tree through rho(s), the
Perron-Frobenius eigenvalue of a matrix m(s) indexed by colours, and state
that Volkov's own 5x+1 construction is a special case of that model. This
script builds the matrix that the TRUE arithmetic tree induces, with the
node's residue as its type, and asks whether its Perron root reproduces
alpha_- or moves away from it.

OUTCOME OF THIS SCRIPT, recorded rather than deleted (Rule 7/14: a
documented dead end is what stops the next session repeating it).

Its central check, `check_type_map`, asks whether the children's residues
mod q^K are a function of the parent's residue mod q^K. That FAILED for
every q and every K tried, so this script never produced a usable matrix,
and `build_matrix` below is unreachable in practice.

The failure is structural, not a bug: w = (2^a u - 1)/q divides by q, so a
congruence on u modulo q^K only pins the child modulo q^(K-1). Demanding
closure at the SAME modulus was the wrong question. The right one is graded
closure, u mod q^(K+1) determining the child mod q^K, and that does hold at
every level; see `type_closure.py`, which supersedes this file and carries
the result.

Method, exactly:
  type of a node   = its residue mod q^K (K given; K=2,3 both run)
  children of u    = w_j = (2^(a0 + j d) u - 1)/q,  j >= 0,  a0 = A0(u mod q)
  displacement     = (a0 + j d) - log2(q)     [in log2 of the value]
  matrix entry     m(s)[r][r'] = sum over children of type r' of 2^(-s*disp)
  growth exponent  = the s where the Perron root of m(s) equals 1

Only FERTILE types carry the recursion (a sterile node is counted but has no
children), so the matrix is restricted to fertile residues, which is what
makes it a genuine multitype object rather than the scalar average.

The type map is VERIFIED to be well defined at each K before any matrix is
built: if the children's residues mod q^K were not a function of the
parent's residue mod q^K, the whole construction would be invalid, and the
script says so and stops rather than reporting a number.
"""
import sys
import numpy as np


def ord_mod(a, m):
    d, x = 1, a % m
    while x != 1:
        x = (x * a) % m
        d += 1
        if d > m:
            raise RuntimeError("no order")
    return d


def a0_of(r, q, d):
    """unique a in 1..d with 2^a r == 1 (mod q); None if r is sterile."""
    for a in range(1, d + 1):
        if (pow(2, a, q) * r) % q == 1 % q:
            return a
    return None


def children(u, q, d, nmax):
    """the first nmax children of odd u in the accelerated reverse tree."""
    a0 = a0_of(u % q, q, d)
    if a0 is None:
        return []
    out = []
    for j in range(nmax):
        w = (pow(2, a0 + j * d) * u - 1) // q
        out.append((a0 + j * d, w))
    return out


def check_type_map(q, d, K, ntest=4000):
    """is (exponent, child residue mod q^K) a function of u mod q^K alone?"""
    M = q ** K
    seen = {}
    bad = 0
    u = 1
    tested = 0
    while tested < ntest:
        u += 2
        if u % q == 0:
            continue
        if a0_of(u % q, q, d) is None:
            continue
        sig = tuple((a, w % M) for a, w in children(u, q, d, 2 * q))
        key = u % M
        if key in seen:
            if seen[key] != sig:
                bad += 1
        else:
            seen[key] = sig
        tested += 1
    return bad, len(seen)


def build_matrix(q, d, K, s, jmax):
    """m(s) over fertile residues mod q^K. Entry = sum of 2^(-s*displacement)."""
    M = q ** K
    types = [r for r in range(M) if r % 2 == 1 and a0_of(r % q, q, d) is not None]
    idx = {r: i for i, r in enumerate(types)}
    n = len(types)
    m = np.zeros((n, n))
    for r in types:
        # a representative odd integer with this residue. NOTE: the check
        # in check_type_map FAILED for every q tried, so this is NOT
        # certified and build_matrix below was never validly usable.
        u = r if r % 2 == 1 else r + M
        while u < 3:
            u += M
        for a, w in children(u, q, d, jmax):
            rw = w % M
            if rw % 2 == 0 or a0_of(rw % q, q, d) is None:
                continue  # sterile child: counted, but carries no recursion
            disp = a - np.log2(q)
            m[idx[r], idx[rw]] += 2.0 ** (-s * disp)
    return m, types


def perron(m):
    ev = np.linalg.eigvals(m)
    return float(np.max(ev.real))


def alpha_minus(q):
    lo, hi = 1e-9, 1.0 - 1e-12
    f = lambda a: (q ** (a - 1.0)) - (2.0 ** a - 1.0)
    for _ in range(400):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def solve_root(q, d, K, jmax):
    """the s at which the Perron root of m(s) equals 1."""
    lo, hi = 0.30, 0.99
    flo = perron(build_matrix(q, d, K, lo, jmax)[0]) - 1.0
    fhi = perron(build_matrix(q, d, K, hi, jmax)[0]) - 1.0
    if flo * fhi > 0:
        return None, (flo, fhi)
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        fm = perron(build_matrix(q, d, K, mid, jmax)[0]) - 1.0
        if flo * fm <= 0:
            hi, fhi = mid, fm
        else:
            lo, flo = mid, fm
    return 0.5 * (lo + hi), None


def main():
    qs = [int(x) for x in sys.argv[1:]] or [3, 5, 7]
    for q in qs:
        d = ord_mod(2, q)
        am = alpha_minus(q)
        print(f"\n=== q={q}  d=ord_q(2)={d}  scalar alpha_- = {am:.9f} ===")
        c = ((2 ** d - 1) // q) % q
        print(f"    H-002 constant c = ((2^d-1)/q) mod q = {c}"
              f"   {'(invertible)' if c % q else '(ZERO: Wieferich)'}")
        for K in (2, 3):
            bad, ntypes = check_type_map(q, d, K)
            status = "well defined" if bad == 0 else f"NOT well defined ({bad} clashes)"
            print(f"    K={K}: type map mod q^{K}={q**K} is {status}"
                  f"  ({ntypes} residues seen)")
            if bad:
                print("    -> refusing to build a matrix on an ill-defined type map")
                continue
            for jmax in (4 * q, 8 * q, 16 * q):
                root, info = solve_root(q, d, K, jmax)
                if root is None:
                    print(f"      jmax={jmax:3d}: no sign change, perron-1 at ends = {info}")
                else:
                    m, types = build_matrix(q, d, K, root, jmax)
                    print(f"      jmax={jmax:3d}: matrix {len(types)}x{len(types)}"
                          f"   Perron root at s -> 1 gives s = {root:.9f}"
                          f"   (s - alpha_-) = {root - am:+.9f}")


if __name__ == "__main__":
    main()
