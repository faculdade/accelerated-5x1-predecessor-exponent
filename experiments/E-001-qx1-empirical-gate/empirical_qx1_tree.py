import math, random

# ============================================================
# Arvore reversa real de qx+1: filhos impares de v sao
# w = (2^a v - 1)/q para a>=1 com 2^a v == 1 mod q.
# DFS com pilha (cada no tem 1 pai unico -> sem revisitas se a raiz
# nao pertence a um ciclo). Conta nos impares <= n_max.
# Previsoes (pressao escalar q^{a-1}=2^a-1):
#   q=5: N_v(x) ~ W_v x^0.650919 ; cauda de W_v: expoente 1/0.650919=1.5363
#   q=7: N_v(x) ~ W_v x^0.373501
# ============================================================

def count_tree(q, root, n_max, checkpoints=None):
    """DFS; retorna contagem total e contagens em checkpoints (<=x)."""
    ordq = 1
    x = 2 % q
    while x != 1:
        x = (2 * x) % q
        ordq += 1
    cps = sorted(checkpoints) if checkpoints else []
    counts = [0] * len(cps)
    total = 0
    stack = [root]
    while stack:
        v = stack.pop()
        total += 1
        for i, cp in enumerate(cps):
            if v <= cp:
                counts[i] += 1
        # filhos: menor a valido
        a0 = None
        p = 2 % q
        for a in range(1, ordq + 1):
            if (p * v) % q == 1:
                a0 = a
                break
            p = (p * 2) % q
        if a0 is None:
            continue  # classe sem filhos (r fora de <2>) ou q | v
        a = a0
        while True:
            w = ((1 << a) * v - 1) // q
            if w > n_max:
                break
            if w != root:
                stack.append(w)
            a += ordq
    return total, counts

CYCLES = {5: {1, 3, 13, 33, 83, 17, 27, 43}, 7: {1}}

# ---------- Parte 1: expoente de crescimento por decada ----------
print("=== Parte 1: expoente local de crescimento, q=5 (previsto 0.6509) ===")
cps = [10**e for e in range(3, 11)]
for root in [7, 11, 19]:
    tot, counts = count_tree(5, root, 10**10, cps)
    slopes = []
    for i in range(1, len(cps)):
        if counts[i - 1] > 0:
            slopes.append(math.log10(counts[i] / counts[i - 1]))
    print(f"raiz {root}: N(1e10)={counts[-1]}, slopes/decada = "
          + " ".join(f"{s:.3f}" for s in slopes))

print("\n=== q=7 (previsto 0.3735) ===")
cps7 = [10**e for e in range(4, 13)]
for root in [9, 11, 15]:
    tot, counts = count_tree(7, root, 10**12, cps7)
    slopes = []
    for i in range(1, len(cps7)):
        if counts[i - 1] > 0:
            slopes.append(math.log10(counts[i] / counts[i - 1]))
    print(f"raiz {root}: N(1e12)={counts[-1]}, slopes/decada = "
          + " ".join(f"{s:.3f}" for s in slopes))

# ---------- Parte 2: cauda de W_v via Hill, q=5 ----------
print("\n=== Parte 2: Hill sobre W_v = N_v(vH)/H^0.650919, q=5, H=1e6 ===")
random.seed(123)
ALPHA1 = 0.650918639898
H = 10**6
sample = []
tried = 0
while len(sample) < 600 and tried < 20000:
    tried += 1
    v = random.randrange(1001, 30001, 2)
    if v % 5 == 0 or v in CYCLES[5]:
        continue
    tot, _ = count_tree(5, v, v * H)
    W = tot / (H ** ALPHA1)
    sample.append(W)
sample.sort(reverse=True)
n = len(sample)
print(f"n={n}, W: min={sample[-1]:.4f} mediana={sample[n//2]:.4f} max={sample[0]:.4f}")
for frac in [0.02, 0.05, 0.10]:
    k = max(5, int(n * frac))
    xk = sample[k]
    hill = k / sum(math.log(sample[i] / xk) for i in range(k))
    print(f"Hill top {frac*100:.0f}% (k={k}): alpha_cauda = {hill:.3f}   (previsto 1.536)")
