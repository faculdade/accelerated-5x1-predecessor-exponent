#!/usr/bin/env python3
"""
E-003, part 2: WHERE the residue type space fails to close, and what the
matrix gives once the failure is quantified.

Part 1 found that "child residue mod q^K is a function of parent residue
mod q^K" is false for every K tried. This script separates three questions
that the blunt check conflated:

  (Q1) graded closure: does u mod q^(K+1) determine each child's residue
       mod q^K? That is the shape the division by q predicts, since
       w = (2^a u - 1)/q loses exactly one factor of q of precision.
  (Q2) fertility closure: does u mod q^2 determine the FERTILITY PATTERN
       of the children, i.e. which sibling positions j are sterile? That
       is all the growth rate needs, and it is weaker than Q1.
  (Q3) equidistribution: given the parent's type, are the children's types
       uniform over the admissible classes? If they are, the annealed
       scalar treatment is exact for the exponent and no finite matrix can
       beat it; if they are not, a matrix could shift the exponent.

Then, under the fertility typing of Q2, it builds m(s) with the child-type
law measured empirically (Q3) rather than assumed, and reports the Perron
root against the scalar alpha_-.
"""
import os
import sys
import numpy as np
from collections import defaultdict, Counter


def ord_mod(a, m):
    d, x = 1, a % m
    while x != 1:
        x = (x * a) % m; d += 1
    return d


def a0_of(r, q, d):
    for a in range(1, d + 1):
        if (pow(2, a, q) * r) % q == 1 % q:
            return a
    return None


def kids(u, q, d, n):
    a0 = a0_of(u % q, q, d)
    if a0 is None:
        return []
    return [(a0 + j * d, (pow(2, a0 + j * d) * u - 1) // q) for j in range(n)]


def alpha_minus(q):
    if q == 3:
        return 1.0
    lo, hi = 1e-9, 1.0 - 1e-12
    f = lambda a: (q ** (a - 1.0)) - (2.0 ** a - 1.0)
    for _ in range(400):
        m = 0.5 * (lo + hi)
        if f(m) > 0: lo = m
        else: hi = m
    return 0.5 * (lo + hi)


def sample_us(q, d, n):
    out, u = [], 1
    while len(out) < n:
        u += 2
        if u % q and a0_of(u % q, q, d) is not None:
            out.append(u)
    return out


def q1_graded(q, d, K, us):
    """does u mod q^(K+1) determine each child's residue mod q^K ?"""
    hi, lo = q ** (K + 1), q ** K
    seen, bad = {}, 0
    for u in us:
        sig = tuple((a, w % lo) for a, w in kids(u, q, d, 2 * q))
        k = u % hi
        if k in seen:
            bad += (seen[k] != sig)
        else:
            seen[k] = sig
    return bad, len(seen)


def q2_fertility(q, d, us):
    """does u mod q^2 determine WHICH sibling positions are sterile?"""
    seen, bad = {}, 0
    for u in us:
        pat = tuple(a0_of(w % q, q, d) is not None for _, w in kids(u, q, d, 3 * q))
        k = u % (q * q)
        if k in seen:
            bad += (seen[k] != pat)
        else:
            seen[k] = pat
    return bad, len(seen)


def q3_equidist(q, d, us):
    """given parent type mod q^2, how are the children's types mod q^2 spread?"""
    tab = defaultdict(Counter)
    for u in us:
        for _, w in kids(u, q, d, 2 * q):
            if a0_of(w % q, q, d) is not None:
                tab[u % (q * q)][w % (q * q)] += 1
    devs = []
    for pt, cnt in tab.items():
        tot = sum(cnt.values())
        k = len(cnt)
        if k < 2:
            continue
        exp = tot / k
        devs.append(max(abs(v - exp) for v in cnt.values()) / exp)
    return (len(tab), (sum(devs) / len(devs) if devs else 0.0),
            (max(devs) if devs else 0.0))


def matrix_root(q, d, us, jmax):
    """m(s) over fertility types mod q^2, child-type law MEASURED, not assumed."""
    T = defaultdict(lambda: defaultdict(float))
    W = Counter()
    for u in us:
        pt = u % (q * q)
        W[pt] += 1
        for a, w in kids(u, q, d, jmax):
            if a0_of(w % q, q, d) is None:
                continue
            T[pt][(w % (q * q), a)] += 1.0
    types = sorted(W)
    idx = {t: i for i, t in enumerate(types)}
    lg = np.log2(q)

    def perron(s):
        M = np.zeros((len(types), len(types)))
        for pt in types:
            for (ct, a), c in T[pt].items():
                if ct in idx:
                    M[idx[pt], idx[ct]] += (c / W[pt]) * 2.0 ** (-s * (a - lg))
        return float(np.max(np.linalg.eigvals(M).real))

    # rho(s) is convex here (negative displacements exist), so rho(s)=1 has
    # TWO roots and the interval's endpoints can share a sign while both roots
    # sit inside. Testing only the endpoints reported "no root" for q=11,23,31
    # when in fact rho dips to 0.73, 0.48 and 0.39 in the interior. Scan first,
    # bracket on [lo, argmin], which is where the SMALLER root lives.
    import numpy as _np
    _grid = _np.linspace(1e-4, 0.999, 121)
    _vals = [perron(float(s)) for s in _grid]
    _imin = int(_np.argmin(_vals))
    lo, hi = 1e-4, float(_grid[_imin])
    flo, fhi = perron(lo) - 1, perron(hi) - 1
    # RANK is the diagnostic that matters more than the root. If the matrix
    # collapses to rank d on a type space of size q(q-1), it is a re-encoding
    # of the scalar pressure and its Perron root carries no new information.
    Mr = _np.zeros((len(types), len(types)))
    for pt in types:
        for (ct, a), c in T[pt].items():
            if ct in idx:
                Mr[idx[pt], idx[ct]] += (c / W[pt]) * 2.0 ** (-0.5 * (a - lg))
    rank = int(_np.linalg.matrix_rank(Mr, tol=1e-10))
    if flo * fhi > 0:
        return None, len(types), (flo, fhi), rank
    for _ in range(90):
        mid = 0.5 * (lo + hi); fm = perron(mid) - 1
        if flo * fm <= 0: hi, fhi = mid, fm
        else: lo, flo = mid, fm
    return 0.5 * (lo + hi), len(types), None, rank


def _report_q(q):
    """Todo o relatorio de UM q. Em nivel de modulo para o multiprocessing
    conseguir picklá-la; cada q e independente dos outros."""
    d = ord_mod(2, q)
    am = alpha_minus(q)
    us = sample_us(q, d, 20000)
    print(f"\n=== q={q}  d={d}  scalar alpha_- = {am:.9f}  ({len(us)} roots sampled) ===")
    print("  Q1 graded closure: u mod q^(K+1) -> child mod q^K")
    for K in (1, 2, 3):
        bad, n = q1_graded(q, d, K, us)
        print(f"     K={K}: {'CLOSES' if bad == 0 else f'FAILS ({bad} clashes)'}"
              f"   ({n} parent classes mod {q**(K+1)})")
    bad, n = q2_fertility(q, d, us)
    print(f"  Q2 fertility pattern from u mod q^2={q*q}: "
          f"{'CLOSES' if bad == 0 else f'FAILS ({bad} clashes)'}  ({n} classes)")
    nt, mdev, xdev = q3_equidist(q, d, us)
    print(f"  Q3 child-type spread given parent type: {nt} parent types, "
          f"mean max-deviation {mdev:.4f}, worst {xdev:.4f}  (0 = perfectly uniform)")
    # The sample is a perfect residue system mod q^3 IF AND ONLY IF n is a
    # multiple of d*q^2. On such an n the matrix collapses to rank d and its
    # Perron root is the scalar pressure identically, so agreement is a
    # tautology, not a check. Every n originally used here was such a
    # multiple for q=5 (d*q^2 = 100), which is why the first version of this
    # experiment reported twelve-decimal "confirmation". Both cases print.
    period = d * q * q
    n_bal = period * max(1, 20000 // period)
    for label, n in ((f"BALANCED   n={n_bal} (multiple of d*q^2={period})", n_bal),
                     (f"unbalanced n={n_bal + 37}", n_bal + 37)):
        us2 = sample_us(q, d, n)
        root, nty, info, rank = matrix_root(q, d, us2, 8 * q)
        tell = f"rank={rank} of {nty} types (d={d})"
        if root is None:
            print(f"  {label}: NO ROOT, ends {info}   {tell}")
        else:
            print(f"  {label}: s = {root:.12f}  (s - alpha_-) = {root - am:+.3e}   {tell}")


def _capture(q):
    import contextlib
    import io as _io
    buf = _io.StringIO()
    with contextlib.redirect_stdout(buf):
        _report_q(q)
    return buf.getvalue()


def main():
    """Paralelo sobre q, uma unidade por processo.

    Rodada 15 (R15-01): este experimento carregava a TERCEIRA razao errada
    para rodar serial. A rodada 9 mediu 4.4 s e concluiu estar uma ordem de
    grandeza abaixo da linha de um minuto, mas mediu o script SEM
    ARGUMENTOS, enquanto o bloco de evidencia declara
    `type_closure.py 5 7 11 13 23 31`, que leva 288 s, cinco vezes ACIMA
    da linha. Foi medida a invocacao errada.

    A forma e uma varredura de valores independentes, que este
    projeto classifica como sempre paralela. A saida e
    escrita na ordem dos argumentos, entao o log reproduz identico.
    """
    qs = [int(x) for x in sys.argv[1:]] or [3, 5, 7]
    if len(qs) < 2:
        for q in qs:
            _report_q(q)
        return
    import multiprocessing
    with multiprocessing.Pool(min(len(qs), os.cpu_count() or 1)) as pool:
        for chunk in pool.map(_capture, qs, chunksize=1):
            sys.stdout.write(chunk)


if __name__ == "__main__":
    main()
