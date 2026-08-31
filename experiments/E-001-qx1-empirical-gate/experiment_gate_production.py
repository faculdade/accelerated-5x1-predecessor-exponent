#!/usr/bin/env python3
"""
H-001 - portao estatistico de producao: a disputa Kontorovich-Lagarias
(2009) vs. Volkov sobre o expoente de contagem da arvore reversa de
5x+1 (0,650919 vs. 0,678, Delta=0,027). Especificacao exata e regra de
admissibilidade derivadas pelo Fable (ver H-109/H-112/H-001), pilotada
em pilot_gate_5x1.py/pilot2_gate_5x1.py neste mesmo diretorio.

Regra de admissibilidade q=5 (nao ha condicao de paridade, so
integralidade mod 5): no impar u, filhos w=(2^a*u-1)/5 validos sse
2^a*u == 1 (mod 5) <=> a == A0[u mod 5] (mod 4), com
A0={1:4, 2:3, 3:1, 4:2}. u==0 mod 5 nao tem filhos (no esteril, conta
mas nao ramifica). O filho e automaticamente impar quando inteiro.

Correcao anti-truncamento (mesma classe do bug historico de E-018):
o ramo a=1 (u=3 mod 5) ENCOLHE (w~0.4u), entao nos <= checkpoint podem
ter ancestrais ACIMA do checkpoint - o teto de BUSCA (search_bound)
deve ficar >=5 decadas acima do maior checkpoint usado na estimativa,
senao o slope sai enviesado para baixo (~-0.012 na fronteira buffer
3->4, medido no piloto). Aqui: search_bound=1e13, janela de checkpoints
1e5..1e8 (buffer=5 decadas).

Todo vies de truncamento conhecido e NEGATIVO (subconta, empurra o
slope para baixo) - por isso o teste e' enquadrado como unilateral:
o numero medido + qualquer correcao residual so se afasta de 0,678
(Volkov) e se aproxima de 0,650919 (Kontorovich-Lagarias/pressao).

Reproduzir: python3 experiment_gate_production.py
"""
import math
import random
import time
from bisect import bisect_right

A0 = {1: 4, 2: 3, 3: 1, 4: 2}
CYCLE5 = {1, 3, 13, 33, 83, 17, 27, 43}

N_ROOTS = 300
SEARCH_BOUND = 10 ** 13
CPS = [10 ** e for e in range(4, 9)]  # checkpoints 1e4..1e8
WINDOW_LO_IDX, WINDOW_HI_IDX = 1, 4  # 1e5 (idx1) .. 1e8 (idx4) = 3 decadas


def decade_counts(root, cps, search_bound):
    """DFS por pilha, teto de BUSCA separado do teto de CONTAGEM
    (checkpoints) - evita o bug de truncamento do braco a=1 encolhendo."""
    bins = [0] * (len(cps) + 1)
    stack = [root]
    while stack:
        u = stack.pop()
        bins[bisect_right(cps, u)] += 1
        r = u % 5
        if r == 0:
            continue
        a = A0[r]
        while True:
            w = ((u << a) - 1) // 5
            if w > search_bound:
                break
            if w != root:
                stack.append(w)
            a += 4
    out, acc = [], 0
    for b in bins[:-1]:
        acc += b
        out.append(acc)
    return out


def sample_roots(rng, n):
    # raizes bem abaixo do checkpoint inferior da janela (1e5) - uma
    # raiz proxima ou acima da janela degenera a contagem (bug corrigido
    # nesta versao: a amostragem original ia ate 200000, > 1e5)
    roots = []
    seen = set()
    while len(roots) < n:
        v = rng.randrange(101, 9999, 2)
        if v % 5 and v not in CYCLE5 and v not in seen:
            seen.add(v)
            roots.append(v)
    return roots


def bootstrap_ci(vals, n_boot=4000, seed=1):
    rng = random.Random(seed)
    n = len(vals)
    boots = []
    for _ in range(n_boot):
        samp = [vals[rng.randrange(n)] for _ in range(n)]
        boots.append(sum(samp) / n)
    boots.sort()
    return boots[int(0.025 * n_boot)], boots[int(0.975 * n_boot) - 1]


def main():
    print(f"=== Portao estatistico de producao (H-001): {N_ROOTS} raizes, "
          f"search_bound={SEARCH_BOUND:.0e}, janela 1e{4+WINDOW_LO_IDX}..1e{4+WINDOW_HI_IDX} ===")
    print("Alvos: Kontorovich-Lagarias/pressao=0.650919 | Volkov=0.678\n")
    rng = random.Random(2026)
    roots = sample_roots(rng, N_ROOTS)

    t0 = time.time()
    slopes = []
    for i, v in enumerate(roots):
        counts = decade_counts(v, CPS, SEARCH_BOUND)
        if counts[WINDOW_LO_IDX] > 0:
            s = math.log10(counts[WINDOW_HI_IDX] / counts[WINDOW_LO_IDX]) / (WINDOW_HI_IDX - WINDOW_LO_IDX)
            slopes.append(s)
        if (i + 1) % 50 == 0:
            print(f"  [{time.time()-t0:6.1f}s] {i+1}/{N_ROOTS} raizes processadas", flush=True)

    n = len(slopes)
    m = sum(slopes) / n
    sd = (sum((s - m) ** 2 for s in slopes) / (n - 1)) ** 0.5
    se = sd / math.sqrt(n)
    lo, hi = bootstrap_ci(slopes, seed=7)

    print(f"\n=== Resultado (n={n}, tempo={time.time()-t0:.1f}s) ===")
    print(f"slope medio = {m:.5f}  sd={sd:.5f}  SE={se:.5f}")
    print(f"IC 95% (bootstrap por raiz) = [{lo:.5f}, {hi:.5f}]")
    print(f"distancia a 0.678 (Volkov) em SEs: {(0.678-m)/se:.2f}")
    print(f"distancia a 0.650919 (KL/pressao) em SEs: {(0.650919-m)/se:.2f}")

    print("\n=== Painel de convergencia (buffer, Richardson) - 5 raizes ===")
    for v in roots[:5]:
        row = []
        for sb_exp in (9, 10, 11, 12, 13):
            counts = decade_counts(v, CPS, 10 ** sb_exp)
            if counts[WINDOW_LO_IDX] > 0:
                s = math.log10(counts[WINDOW_HI_IDX] / counts[WINDOW_LO_IDX]) / (WINDOW_HI_IDX - WINDOW_LO_IDX)
                row.append(f"1e{sb_exp}:{s:.4f}")
        print(f"  root {v:6d}: " + " ".join(row))


if __name__ == "__main__":
    main()


# Migration note: references `E-018`, `H-109`, `H-112`, out of this repo's scope, left unrenumbered.
