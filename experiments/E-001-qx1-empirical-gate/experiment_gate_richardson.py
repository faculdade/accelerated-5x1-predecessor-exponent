#!/usr/bin/env python3
"""
H-001 - versao com rastreio de path-max, proposta pelo Fable apos a
producao inicial (n=300) dar um IC bilateral que excluia os DOIS
valores candidatos (0.678 e 0.650919) - modo de falha previsto por ele
mesmo: o IC bruto e' um IC do ESTIMADOR enviesado-para-baixo, nao de
gamma, e o vies residual (~0.01) e maior que o SE (~0.003) em n=300.

Em vez de rodar o DFS 5x por raiz (uma por search_bound, como no
painel de convergencia anterior), rastreamos path_max = maior valor de
no visitado ao longo do caminho raiz->no (path_max(filho) =
max(path_max(pai), filho)). Um no so teria sido gerado num run com
search_bound=10^B se TODO o caminho ate ele (raiz inclusive) ficou
<=10^B - ou seja, path_max(no)<=10^B. Isso da, de uma unica passada em
search_bound=1e13, as contagens em TODOS os buffers 9..13
simultaneamente (binning 2D: checkpoint de valor x buffer de path_max,
depois soma cumulativa nas duas dimensoes) - custo computacional igual
a UM DFS grande, nao 5.

Extrapolacao de Richardson (Aitken Delta^2) na curva MEDIA entre raizes
(nunca por raiz - series por raiz sao granulosas/nao-monotonas e o
Aitken diverge; a media pooled decai quase geometricamente, confirmado
no painel anterior) estima o valor no limite buffer->infinito.

Reproduzir: python3 experiment_gate_richardson.py
"""
import math
import multiprocessing
import os
import random
import time
from bisect import bisect_left, bisect_right

A0 = {1: 4, 2: 3, 3: 1, 4: 2}
CYCLE5 = {1, 3, 13, 33, 83, 17, 27, 43}

N_ROOTS = 300
SEARCH_BOUND = 10 ** 13
CPS = [10 ** e for e in range(4, 9)]  # checkpoints de VALOR: 1e4..1e8
BUFFERS = [9, 10, 11, 12, 13]  # buffers (expoente de path_max) testados
WINDOW_LO_IDX, WINDOW_HI_IDX = 1, 4  # janela 1e5..1e8 (3 decadas)


def decade_counts_2d_v2(root, cps, buffers, search_bound):
    """DFS com pilha rastreando path_max (path_max(filho)=max(path_max(pai),
    filho)); acumula contagem bruta por (checkpoint, buffer_exato de
    path_max), depois faz prefix-sum 2D (checkpoint cumulativo, buffer
    cumulativo) para dar N(valor<=cps[ci], path_max<=10^buffers[bi])."""
    n_cp, n_buf = len(cps), len(buffers)
    # limiares INTEIROS, como o enumerador em C (tree_counts.c:215 compara
    # bufthr[bi] < pmax). A versao anterior fazia
    # bisect_right(buffers, log10(pmax) + 1e-9), e o epsilon empurrava os
    # valores logo abaixo de uma potencia de dez para o buffer seguinte:
    # para pmax = 999999999, log10 + 1e-9 passa de 9, entao o no ia para o
    # buffer 1e10 em vez do 1e9. Achado na rodada 5 (R5-15).
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
    # prefix-sum: cumulativo em checkpoint (dim 0) e em buffer (dim 1)
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
    return raw  # raw[ci][bi] = N(valor<=cps[ci], path_max<=10^buffers[bi])


def sample_roots(rng, n):
    roots = []
    seen = set()
    while len(roots) < n:
        v = rng.randrange(101, 9999, 2)
        if v % 5 and v not in CYCLE5 and v not in seen:
            seen.add(v)
            roots.append(v)
    return roots


def aitken(s):
    """Aitken Delta^2 nos ultimos 3 pontos da sequencia s (lista)."""
    if len(s) < 3:
        return None
    s1, s2, s3 = s[-3], s[-2], s[-1]
    denom = (s3 - s2) - (s2 - s1)
    if abs(denom) < 1e-12:
        return None
    return s3 - (s3 - s2) ** 2 / denom


def bootstrap_richardson(all_slopes_by_root, n_boot=2000, seed=3):
    """Bootstrap: reamostra raizes, recalcula curva media por buffer, aplica Aitken."""
    rng = random.Random(seed)
    n = len(all_slopes_by_root)
    n_buf = len(all_slopes_by_root[0])
    boots = []
    for _ in range(n_boot):
        idxs = [rng.randrange(n) for _ in range(n)]
        means = []
        for bi in range(n_buf):
            means.append(sum(all_slopes_by_root[i][bi] for i in idxs) / n)
        a = aitken(means)
        if a is not None:
            boots.append(a)
    boots.sort()
    return boots



def _slopes_for_root(v):
    """One root's work. Module level so multiprocessing can pickle it.
    Independent of every other root, which is why Rule 9c calls this shape
    'independent trials, always parallel'."""
    cum = decade_counts_2d_v2(v, CPS, BUFFERS, SEARCH_BOUND)
    out = []
    for bi in range(len(BUFFERS)):
        lo = cum[WINDOW_LO_IDX][bi]
        hi = cum[WINDOW_HI_IDX][bi]
        if lo > 0 and hi > 0:
            out.append(math.log10(hi / lo) / (WINDOW_HI_IDX - WINDOW_LO_IDX))
        else:
            out.append(float("nan"))
    return out


def main():
    print(f"=== Portao com Richardson (path-max, uma passada): {N_ROOTS} raizes, "
          f"search_bound={SEARCH_BOUND:.0e} ===")
    print(f"Buffers testados: {BUFFERS}  Janela: 1e{4+WINDOW_LO_IDX}..1e{4+WINDOW_HI_IDX}\n")
    rng = random.Random(2026)
    roots = sample_roots(rng, N_ROOTS)

    t0 = time.time()
    # Parallel over roots. Pool.map preserves input order, so the resampling
    # in bootstrap_richardson sees the same list it saw when this ran serially
    # and the output is unchanged. Rule 9c: independent trials, always parallel.
    n_workers = min(os.cpu_count() or 1, N_ROOTS)
    with multiprocessing.Pool(n_workers) as pool:
        all_slopes_by_root = pool.map(_slopes_for_root, roots, chunksize=4)
    print(f"  [{time.time()-t0:6.1f}s] {len(all_slopes_by_root)}/{N_ROOTS} raizes "
          f"em {n_workers} processos", flush=True)

    print(f"\n=== Curva media por buffer (n={len(all_slopes_by_root)}, tempo={time.time()-t0:.1f}s) ===")
    n_buf = len(BUFFERS)
    mean_by_buffer = []
    for bi in range(n_buf):
        vals = [row[bi] for row in all_slopes_by_root if not math.isnan(row[bi])]
        m = sum(vals) / len(vals)
        mean_by_buffer.append(m)
        print(f"  buffer=1e{BUFFERS[bi]}: media={m:.5f}  n={len(vals)}")

    print("\nincrementos:", [f"{mean_by_buffer[i+1]-mean_by_buffer[i]:.5f}" for i in range(n_buf-1)])
    a_last3 = aitken(mean_by_buffer)
    a_prev3 = aitken(mean_by_buffer[:-1])
    print(f"Aitken (ultimos 3 pontos, B11-13): correcao->  s_inf = {a_last3:.5f}" if a_last3 else "Aitken B11-13: indefinido")
    print(f"Aitken (B10-12): s_inf = {a_prev3:.5f}" if a_prev3 else "Aitken B10-12: indefinido")

    print("\n=== Bootstrap do valor extrapolado (Aitken B11-13, por raiz reamostrada) ===")
    boots = bootstrap_richardson(all_slopes_by_root, seed=3)
    lo_ci = boots[int(0.025 * len(boots))]
    hi_ci = boots[int(0.975 * len(boots)) - 1]
    point = sum(boots) / len(boots)
    print(f"n_boot_validos={len(boots)}  media={point:.5f}  IC95%=[{lo_ci:.5f}, {hi_ci:.5f}]")
    print(f"distancia de 0.678 (Volkov) ao limite superior do IC: {0.678 - hi_ci:.4f}")
    print(f"distancia de 0.650919 (KL/pressao) ao centro: {0.650919 - point:.4f}")

    print("\n=== Painel de slopes por decada dentro da janela (buffer=1e13, teste de pre-assintotica) ===")
    print("(se o slope por decada ainda cresce com k, o residuo ate 0.6509 e' pre-assintotico, nao truncamento)")
    # recomputa por decada usando o maior buffer diretamente (reusa cum do ultimo root calculado nao serve;
    # roda um pequeno subconjunto de raizes de novo so para o diagnostico por decada)
    sub_roots = roots[:40]
    per_decade = [[] for _ in range(len(CPS) - 1)]
    for v in sub_roots:
        cum = decade_counts_2d_v2(v, CPS, BUFFERS, SEARCH_BOUND)
        bi_max = len(BUFFERS) - 1
        for ci in range(len(CPS) - 1):
            lo = cum[ci][bi_max]
            hi = cum[ci + 1][bi_max]
            if lo > 0 and hi > 0:
                per_decade[ci].append(math.log10(hi / lo))
    for ci in range(len(CPS) - 1):
        vals = per_decade[ci]
        if vals:
            m = sum(vals) / len(vals)
            print(f"  decada 1e{4+ci}->1e{5+ci}: slope_medio={m:.4f}  n={len(vals)}")


if __name__ == "__main__":
    main()
