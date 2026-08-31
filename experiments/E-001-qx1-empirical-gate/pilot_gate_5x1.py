"""Piloto do portao estatistico H-001: calibrar ruido/vies/custo do
slope de contagem por decada na arvore reversa REAL de 5x+1.

Regra exata (q=5): no impar u, filhos w=(2^a*u - 1)/5, validos sse
2^a*u == 1 (mod 5)  <=>  a == a0(u mod 5) (mod 4), com
a0: {1:4, 2:3, 3:1, 4:2}; u==0 mod 5 nao tem filhos.
Filho e automaticamente impar (5w = 2^a*u - 1 impar => w impar).

Correcao anti-truncamento (bug do E-018): limite de BUSCA (search_bound)
separado do limite de CONTAGEM (checkpoints). Ramos a=1 (u=3 mod 5)
ENCOLHEM (w ~ 0.4u), entao nos <= checkpoint podem ter ancestrais acima
do checkpoint; medimos a sensibilidade variando search_bound.
"""
import math, time, random
from bisect import bisect_right

A0 = {1: 4, 2: 3, 3: 1, 4: 2}

def decade_counts(root, cps, search_bound):
    """Retorna N(cp) cumulativo para cada checkpoint, DFS ate search_bound."""
    bins = [0] * (len(cps) + 1)
    stack = [root]
    nodes = 0
    while stack:
        u = stack.pop()
        nodes += 1
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
    # cumulativo
    out, acc = [], 0
    for b in bins[:-1]:
        acc += b
        out.append(acc)
    return out, nodes

CYCLE5 = {1, 3, 13, 33, 83, 17, 27, 43}

cps = [10**e for e in range(2, 10)]  # 1e2..1e9
random.seed(42)
roots = []
while len(roots) < 8:
    v = random.randrange(7, 500, 2)
    if v % 5 and v not in CYCLE5 and v not in roots:
        roots.append(v)

print("== sensibilidade ao search_bound (vies de truncamento) ==")
print("root | bound | N(1e6) N(1e8) N(1e9) | nos_total | t(s)")
for root in roots[:3]:
    for sb_exp in (9, 10, 11):
        t0 = time.time()
        counts, nodes = decade_counts(root, cps, 10**sb_exp)
        dt = time.time() - t0
        print(f"{root:4d} | 1e{sb_exp} | {counts[4]:7d} {counts[6]:8d} "
              f"{counts[7]:9d} | {nodes:9d} | {dt:.1f}")

print()
print("== slopes por decada (search_bound=1e11, contagem ate 1e9) ==")
print("previsto KL/pressao: 0.650919 ; Volkov: 0.678")
slopes_all = []
for root in roots:
    t0 = time.time()
    counts, nodes = decade_counts(root, cps, 10**11)
    dt = time.time() - t0
    sl = [math.log10(counts[i] / counts[i - 1])
          for i in range(1, len(cps)) if counts[i - 1] > 0]
    # janela estavel: decadas 1e5..1e9 (4 ultimos slopes)
    stable = sl[-4:]
    slopes_all.extend(stable)
    print(f"root {root:4d}: N(1e9)={counts[-1]:8d} slopes="
          + " ".join(f"{s:.4f}" for s in sl) + f"  t={dt:.1f}s")

m = sum(slopes_all) / len(slopes_all)
sd = (sum((s - m) ** 2 for s in slopes_all) / (len(slopes_all) - 1)) ** 0.5
print(f"\njanela estavel (1e5..1e9): media={m:.4f} sd={sd:.4f} "
      f"n={len(slopes_all)} SE={sd/len(slopes_all)**0.5:.4f}")


# Migration note: references `E-018`, out of this repo's scope, left unrenumbered.
