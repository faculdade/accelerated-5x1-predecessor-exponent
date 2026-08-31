import numpy as np
from math import log, gcd

# ============================================================
# Pressao multitipo para arvores reversas de qx+1 (q impar).
# Recursao: G(v) = sum_a q*2^{-a} G(w_a), w_a = (2^a v - 1)/q,
# a valido sse 2^a v == 1 mod q, a>=1.
# Tipos = residuos r mod q^k "vivos" (gcd(r,q)=1 E existe a valido,
# i.e. r mod q pertence ao subgrupo <2> de (Z/q)^*).
# M(alpha)_{r,s} = sum_{a valido} (q 2^{-a})^alpha * P_d[filho = s vivo],
# d = proximo digito q-adico, uniforme em {0..q-1}.
# TEOREMA-CANDIDATO (somas de coluna constantes):
#   rho(M_q,k(alpha)) = c_q(alpha) := (1/q) q^alpha 2^{-alpha}/(1-2^{-alpha})
#                     = q^{alpha-1}/(2^alpha - 1)   para todo k.
# ============================================================

def subgroup_of_2(q):
    s, x = set(), 1
    while True:
        x = (x * 2) % q
        if x in s:
            break
        s.add(x)
    return s

def build_matrix(q, k, alpha, A=400):
    mod = q ** k
    H2 = subgroup_of_2(q)
    alive = [r for r in range(mod) if gcd(r, q) == 1 and (pow(r, -1, q) in H2)]
    idx = {r: i for i, r in enumerate(alive)}
    n = len(alive)
    M = np.zeros((n, n))
    aliveset = set(alive)
    for r in alive:
        for a in range(1, A + 1):
            if (pow(2, a, q) * r) % q != 1:
                continue
            w_weight = (q * 2.0 ** (-a)) ** alpha
            for d in range(q):
                v = r + d * mod
                num = (2 ** a) * v - 1
                assert num % q == 0
                w = (num // q) % mod
                if w in aliveset:
                    M[idx[r], idx[w]] += w_weight / q
    return M

def c_scalar(q, alpha):
    return q ** (alpha - 1) / (2 ** alpha - 1)

print("=== verificacao: rho(M_q,k(alpha)) vs forma fechada escalar ===")
for q in [3, 5, 7]:
    kmax = {3: 4, 5: 3, 7: 3}[q]
    for k in range(1, kmax + 1):
        for alpha in [0.5, 1.0, 1.5, 2.0]:
            M = build_matrix(q, k, alpha)
            rho = max(abs(np.linalg.eigvals(M)))
            c = c_scalar(q, alpha)
            # tambem checa somas de coluna
            cs = M.sum(axis=0)
            print(f"q={q} k={k} alpha={alpha:.1f}: rho={rho:.12f} "
                  f"c={c:.12f} diff={abs(rho-c):.2e} "
                  f"colsum[min,max]=[{cs.min():.10f},{cs.max():.10f}]")
    print()

# ============================================================
# raizes de q^(alpha-1) = 2^alpha - 1
# ============================================================
from scipy.optimize import brentq

def f(alpha, q):
    return (alpha - 1) * log(q) - log(2 ** alpha - 1)

print("=== raizes da equacao de pressao q^(a-1) = 2^a - 1 ===")
print(f"{'q':>3} {'raiz menor a1':>16} {'raiz maior a2':>16} {'cauda k=a2/a1':>16}")
for q in [3, 5, 7, 9, 11, 13, 15]:
    # alpha=1 e sempre raiz; achar a outra
    roots = [1.0]
    # procura em (0.01, 0.999) e (1.001, 6)
    for lo, hi in [(0.01, 0.999), (1.001, 6.0)]:
        try:
            if f(lo, q) * f(hi, q) < 0:
                roots.append(brentq(f, lo, hi, args=(q,), xtol=1e-14))
        except ValueError:
            pass
    roots = sorted(roots)
    a1, a2 = roots[0], roots[-1]
    print(f"{q:>3} {a1:>16.12f} {a2:>16.12f} {a2/a1:>16.12f}")
