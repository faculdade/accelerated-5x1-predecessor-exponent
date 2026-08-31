#!/bin/sh
# E-004: reproduz out/volkov_cutoff.log inteiro, cabecalho incluido.
#
# Escrito na rodada 14 (R14-06). Ate aqui o bloco de evidencia trazia
# `./volkov_count <m> <T> for m in 1e3..1e6 and T in 4.2e8..1e16`, que e um
# ESQUEMA e nao roda como esta escrito, e o log carregava linhas de
# cabecalho compostas a mao que o programa nao emite. Um revisor nao tinha
# como reproduzir o arquivo.
# Locale fixo: sem isto o awk usa virgula decimal em pt_BR e o log sai
# diferente conforme a maquina do revisor. Achado ao escrever este script
# na rodada 14.
LC_ALL=C
export LC_ALL

cd "$(dirname "$0")" || exit 1
gcc -O2 -o volkov_count volkov_count.c -lm || exit 1

cat <<'HDR'
# E-004, mapa CORRETO da eq (1.1) do Volkov: M(x)=5x+1 se impar, x/2 se par
# Ciclos calculados a partir de 1, 13, 17 (27 valores, o primeiro completado:
# o texto dele imprime 1->6->3->1, mas M(3)=16)

## Reproducao das quatro estimativas impressas (ele: 0.741, 0.715, 0.694, 0.680)
HDR
for m in 1000 10000 100000 1000000; do ./volkov_count "$m" 420000000; done

echo
echo "## Dependencia do limiar, m fixo em 1e6"
for T in 420000000 1000000000 10000000000 1000000000000 10000000000000000; do
  ./volkov_count 1000000 "$T"
done

echo
echo "## Dependencia do limiar em CADA m (acrescentado 2026-08-29, rodada 5)"
echo "## O paper afirmava que \"os mesmos quatro valores carregam o erro do limiar"
echo "## em todo m\". Instanciado (Regra 11c, passo 3), e falso no extremo inferior."
for m in 1000 10000 100000 1000000; do
  a=$(./volkov_count "$m" 420000000            | sed 's/.*beta=//')
  b=$(./volkov_count "$m" 10000000000000000    | sed 's/.*beta=//')
  d=$(awk -v x="$a" -v y="$b" 'BEGIN{printf "%+.5f", y-x}')
  printf 'm=%-7s T=4.2e8 beta=%s   T=1e16 beta=%s   diferenca %s\n' \
         "$m" "$a" "$b" "$d"
done
echo "## O erro e zero no m menor e cresce monotonamente com m, entao ele nao"
echo "## contamina os quatro por igual: acentua a queda."
