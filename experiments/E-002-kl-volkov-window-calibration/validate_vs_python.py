#!/usr/bin/env python3
"""
E-002 validation: the C enumerator must reproduce, entry by entry, the
count matrix of the original E-001 Python enumerator
(`experiment_gate_richardson.decade_counts_2d_v2`).

Any deviation invalidates every deeper run, so this is run before anything
else. Two earlier stages of this project were derailed by exactly this
class of enumeration bug (the E-018 pruning bug and the H-001 root-sampling
bug), which is why the check is a separate file and not a comment.

Run: python3 validate_vs_python.py
"""
import math
import subprocess
import sys
from bisect import bisect_left, bisect_right

HERE = __file__.rsplit("/", 1)[0]

A0 = {1: 4, 2: 3, 3: 1, 4: 2}


def decade_counts_2d_v2(root, cps, buffers, search_bound):
    """Verbatim from E-001 experiment_gate_richardson.py.

    A rodada 13 achou que esta copia NAO era mais literal: ela ainda usava
    `bisect_right(buffers, log10(pmax) + 1e-9)`, o agrupamento que a rodada
    5 (R5-15) consertou no original, enquanto a docstring seguia dizendo
    "verbatim". Como este script e o portao que valida o enumerador em C
    contra o Python, ele vinha comparando o C contra um Python defeituoso.
    Sincronizado com o original.
    """
    n_cp, n_buf = len(cps), len(buffers)
    bufthr = [10 ** b for b in buffers]
    raw = [[0] * n_buf for _ in range(n_cp)]
    stack = [(root, root)]
    while stack:
        u, pmax = stack.pop()
        ci = bisect_right(cps, u)
        if ci < n_cp:
            bi = bisect_left(bufthr, pmax)
            if bi < n_buf:
                raw[ci][bi] += 1
        r = u % 5
        if r == 0:
            continue
        a = A0[r]
        while True:
            w = ((u << a) - 1) // 5
            if w > search_bound:
                break
            if w != root:
                new_pmax = w if w > pmax else pmax
                stack.append((w, new_pmax))
            a += 4
    for bi in range(n_buf):
        acc = 0
        for ci in range(n_cp):
            acc += raw[ci][bi]
            raw[ci][bi] = acc
    for ci in range(n_cp):
        acc = 0
        for bi in range(n_buf):
            acc += raw[ci][bi]
            raw[ci][bi] = acc
    return raw


def main():
    roots = [103, 1237, 4441, 7919, 9973]
    cps = [10 ** e for e in range(4, 9)]
    buffers = [9, 10, 11, 12, 13]
    bound = 10 ** 13

    cmd = [HERE + "/tree_counts", "--q", "5", "--cp", "4", "8", "--buf", "9", "13",
           "--rootlist", ",".join(str(r) for r in roots)]
    out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout

    cmat = {}
    for line in out.splitlines():
        if line.startswith("#"):
            continue
        f = line.split()
        root = int(f[0])
        vals = [int(x) for x in f[2:]]
        cmat[root] = [vals[ci * len(buffers):(ci + 1) * len(buffers)] for ci in range(len(cps))]

    bad = 0
    for root in roots:
        py = decade_counts_2d_v2(root, cps, buffers, bound)
        c = cmat[root]
        same = (py == c)
        print(f"root {root:6d}: python == C ? {same}   N(1e8,1e13)={py[-1][-1]}")
        if not same:
            bad += 1
            for ci in range(len(cps)):
                if py[ci] != c[ci]:
                    print(f"   cp=1e{4+ci}: python={py[ci]}  C={c[ci]}")
    print("\nVALIDATION", "PASSED" if bad == 0 else f"FAILED on {bad} roots")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())


# Migration note: references `E-018`, out of this repo's scope, left unrenumbered.
