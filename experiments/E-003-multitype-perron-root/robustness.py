#!/usr/bin/env python3
"""E-003, part 3: is the matrix root of part 2 stable, and is the q=7 gap real?

Three stress tests on the part-2 result:
  (a) sample size: does the root move as the number of sampled roots grows?
  (b) sibling truncation jmax: the matrix sums only the first jmax siblings,
      so the root must converge as jmax grows; report the tail.
  (c) is the exact-uniformity of the child-type law an artifact of sampling
      only small u? Re-measure starting from a large offset.
"""
import multiprocessing
import os

import numpy as np
from collections import defaultdict, Counter
from type_closure import ord_mod, a0_of, kids, alpha_minus


def sample_from(q, d, n, start):
    out, u = [], start | 1
    while len(out) < n:
        u += 2
        if u % q and a0_of(u % q, q, d) is not None:
            out.append(u)
    return out


def root_of(q, d, us, jmax):
    T = defaultdict(lambda: defaultdict(float)); W = Counter()
    for u in us:
        pt = u % (q * q); W[pt] += 1
        for a, w in kids(u, q, d, jmax):
            if a0_of(w % q, q, d) is not None:
                T[pt][(w % (q * q), a)] += 1.0
    types = sorted(W); idx = {t: i for i, t in enumerate(types)}; lg = np.log2(q)
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
    if flo * fhi > 0: return None
    for _ in range(90):
        mid = 0.5 * (lo + hi); fm = perron(mid) - 1
        if flo * fm <= 0: hi, fhi = mid, fm
        else: lo, flo = mid, fm
    return 0.5 * (lo + hi)


def _task(job):
    """Uma unidade de trabalho. Em nivel de modulo para o multiprocessing
    conseguir picklá-la. sample_from e root_of sao deterministicos (nao ha
    aleatoriedade: a amostra e a caminhada dos u impares admissiveis a
    partir de start), entao cada worker regenera a sua propria amostra e o
    resultado independe de como o trabalho foi repartido."""
    q, kind, arg = job
    d = ord_mod(2, q)
    if kind == "a":
        return job, root_of(q, d, sample_from(q, d, arg, 1), 8 * q)
    if kind == "b":
        return job, root_of(q, d, sample_from(q, d, 100000, 1), arg)
    return job, root_of(q, d, sample_from(q, d, 50000, 10 ** 12), 16 * q)


def main():
    jobs = []
    for q in (5, 7):
        jobs += [(q, "a", n) for n in (2000, 20000, 100000)]
        jobs += [(q, "b", jm) for jm in (2*q, 4*q, 8*q, 16*q, 32*q)]
        jobs += [(q, "c", None)]

    # Regra 9c: varredura de parametros, uma unidade por processo. O bloco
    # (b) e o mais caro, entao chunksize=1 para nao empilhar os caros num
    # worker so.
    n_workers = min(os.cpu_count() or 1, len(jobs))
    with multiprocessing.Pool(n_workers) as pool:
        res = dict(pool.map(_task, jobs, chunksize=1))

    # impressao na ordem serial original, para a saida ser comparavel
    for q in (5, 7):
        am = alpha_minus(q)
        print(f"\n=== q={q}  alpha_- = {am:.12f} ===")
        print("  (a) sample size, jmax=8q")
        for n in (2000, 20000, 100000):
            r = res[(q, "a", n)]
            print(f"      n={n:6d}: s = {r:.12f}   diff = {r-am:+.12f}")
        print("  (b) sibling truncation jmax, n=100000")
        prev = None
        for jm in (2*q, 4*q, 8*q, 16*q, 32*q):
            r = res[(q, "b", jm)]
            step = "" if prev is None else f"   step = {r-prev:+.3e}"
            print(f"      jmax={jm:4d}: s = {r:.12f}   diff = {r-am:+.12f}{step}")
            prev = r
        print("  (c) large-u window, n=50000 starting at 10^12, jmax=16q")
        r = res[(q, "c", None)]
        print(f"      s = {r:.12f}   diff = {r-am:+.12f}")


if __name__ == "__main__":
    main()
