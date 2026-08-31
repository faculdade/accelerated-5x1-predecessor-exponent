"""Confere a afirmacao do paper de que a passada unica com path-max
concorda CONTAGEM POR CONTAGEM com uma implementacao independente que
re-roda a enumeracao uma vez por buffer.

A rodada 4 da critica achou que essa frase estava no Method sem
nenhuma execucao comitada por tras. As duas implementacoes ja existiam
nesta pasta; faltava compara-las.

  decade_counts_2d_v2  (experiment_gate_richardson.py): uma passada,
      rastreia o maximo do caminho e devolve os cinco buffers de uma vez.
  decade_counts        (experiment_gate_production.py): re-roda a
      enumeracao inteira com um limite duro diferente por buffer.

Sao algoritmos diferentes, nao a mesma rotina em duas linguagens, que e
o que o out/validate_vs_python.log do E-002 ja cobria.
"""
import random, sys, time

from experiment_gate_richardson import (decade_counts_2d_v2, sample_roots,
                                        CPS, BUFFERS, SEARCH_BOUND, N_ROOTS)
from experiment_gate_production import decade_counts

N_CHECK = 5

def main():
    rng = random.Random(2026)
    roots = sample_roots(rng, N_ROOTS)
    print(f"{N_CHECK} das {len(roots)} raizes, as mesmas que o painel de "
          f"convergencia do script de producao usa: {roots[:N_CHECK]}")
    print(f"checkpoints {CPS}")
    print(f"buffers 1e{BUFFERS[0]}..1e{BUFFERS[-1]}\n")

    t0 = time.time()
    bad = 0
    for v in roots[:N_CHECK]:
        single = decade_counts_2d_v2(v, CPS, BUFFERS, SEARCH_BOUND)
        for j, b in enumerate(BUFFERS):
            per = decade_counts(v, CPS, 10 ** b)
            # decade_counts_2d_v2 devolve raw[checkpoint][buffer]
            one = [single[i][j] for i in range(len(CPS))]
            ok = list(per) == list(one)
            if not ok:
                bad += 1
                print(f"  DIVERGE raiz {v} buffer 1e{b}")
                print(f"    passada unica : {one}")
                print(f"    por buffer    : {list(per)}")
            else:
                print(f"  raiz {v:6d} buffer 1e{b}: {one}  igual")
    dt = time.time() - t0
    print(f"\n{N_CHECK} raizes x {len(BUFFERS)} buffers x {len(CPS)} "
          f"checkpoints = {N_CHECK*len(BUFFERS)*len(CPS)} contagens, {dt:.1f}s")
    if bad:
        print(f"FALHOU: {bad} par(es) de buffers divergem")
        return 1
    print("VALIDACAO OK: concorda contagem por contagem")
    return 0

if __name__ == "__main__":
    sys.exit(main())
