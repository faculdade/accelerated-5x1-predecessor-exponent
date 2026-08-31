#!/bin/sh
# E-002: the deep runs. Seven datasets, same roots, same checkpoints, same
# buffers, so every decade is comparable across all of them.
#
# Buffer 1e15 rather than 1e17. The stochastic controls have a heavy-tailed
# total progeny (tail index 1/0.650919 = 1.5363), so at 1e17 a single
# unlucky realization out of 300 dominates the wall time and the batch never
# ends. At 1e15 the deepest fully buffered decade is 1e9 -> 1e10, which is
# two decades deeper than the E-001 window and finishes in minutes.
#
# The arithmetic tree alone was also run to 1e17, since it has no such tail;
# see data/q5_arith_b17.txt and the deep table in the README.
cd "$(dirname "$0")" || exit 1
for m in "" "--cyc" "--iid" "--cycq 5.00000" "--cycq 5.05398" "--iidq 5.00000" "--iidq 5.05398"; do
  case "$m" in
    "")               tag=arith ;;
    "--cyc")          tag=cyc ;;
    "--iid")          tag=iid ;;
    "--cycq 5.00000") tag=cycq500 ;;
    "--cycq 5.05398") tag=cycq505 ;;
    "--iidq 5.00000") tag=iidq500 ;;
    "--iidq 5.05398") tag=iidq505 ;;
  esac
  echo "running $tag"
  # shellcheck disable=SC2086
  ./tree_counts --q 5 $m --roots 300 --cp 4 10 --buf 9 15 --out "data/q5_${tag}_b15.txt"
done

# A arvore aritmetica sozinha tambem foi levada a 1e17, porque ela nao tem a
# cauda pesada dos controles. Este comando faltava no repositorio inteiro
# ate a rodada 14 (R14-04): o arquivo q5_arith_b17.txt sustenta a evidencia
# mais profunda do paper (a tabela 0.6506/0.6505, o "never reads above
# 0.650919" e o termo de truncamento <= 0.0004) e nao havia como um revisor
# saber como regera-lo. Leva cerca de 824 s numa maquina ociosa.
echo "running arith b17 (deep, ~17 min)"
./tree_counts --q 5 --roots 300 --cp 4 12 --buf 9 17 --out data/q5_arith_b17.txt

echo DEEP_DONE
