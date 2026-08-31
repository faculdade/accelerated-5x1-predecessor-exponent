"""Piloto 2: (a) convergencia do slope vs buffer de busca;
(b) estimativa central com 60 raizes, janela pre-registrada, bootstrap."""
import math, time, random
from bisect import bisect_right

A0 = {1: 4, 2: 3, 3: 1, 4: 2}
CYCLE5 = {1, 3, 13, 33, 83, 17, 27, 43}

def decade_counts(root, cps, search_bound):
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

# ---------- (a) slope da decada 1e5->1e6 e 1e6->1e7 vs bound ----------
cps = [10**e for e in range(4, 9)]  # 1e4..1e8
print("== slope fixo vs search_bound (root 333, 77, 19) ==")
print("decada 1e5->1e6 (s56) e 1e6->1e7 (s67); buffer = log10(bound)-7")
for root in (333, 77, 19):
    row = []
    for sb in (8, 9, 10, 11, 12):
        c = decade_counts(root, cps, 10**sb)
        s56 = math.log10(c[2] / c[1])
        s67 = math.log10(c[3] / c[2])
        row.append(f"1e{sb}: s56={s56:.4f} s67={s67:.4f}")
    print(f"root {root:4d}: " + " | ".join(row))

# ---------- (b) 60 raizes, bound 1e12, janela 1e5..1e8 ----------
print("\n== 60 raizes impares em [101,999], bound 1e12, janela 1e5..1e8 ==")
random.seed(7)
roots = []
while len(roots) < 60:
    v = random.randrange(101, 1000, 2)
    if v % 5 and v not in CYCLE5 and v not in roots:
        roots.append(v)

cps2 = [10**e for e in range(4, 9)]  # 1e4..1e8; buffer minimo = 4 decadas
t0 = time.time()
per_root_slope = []
for v in roots:
    c = decade_counts(v, cps2, 10**12)
    # slope medio da janela 1e5..1e8 (3 decadas), por regressao = media
    if c[1] > 0:
        s = math.log10(c[4] / c[1]) / 3.0
        per_root_slope.append(s)
dt = time.time() - t0
n = len(per_root_slope)
m = sum(per_root_slope) / n
sd = (sum((s - m) ** 2 for s in per_root_slope) / (n - 1)) ** 0.5
print(f"n={n} raizes, t={dt:.0f}s: slope medio={m:.4f} sd={sd:.4f} "
      f"SE={sd/n**0.5:.4f}")

# bootstrap percentil (raiz como unidade)
B, boots = 4000, []
for _ in range(B):
    samp = [per_root_slope[random.randrange(n)] for _ in range(n)]
    boots.append(sum(samp) / n)
boots.sort()
lo, hi = boots[int(0.025 * B)], boots[int(0.975 * B)]
print(f"IC 95% bootstrap: [{lo:.4f}, {hi:.4f}]")
print(f"alvos: KL/pressao 0.650919 | Volkov 0.678")
