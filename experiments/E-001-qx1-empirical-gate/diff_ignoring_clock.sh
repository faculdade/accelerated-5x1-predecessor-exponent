#!/bin/bash
# Compara um log recem-gerado com o comitado, neutralizando SO os valores de
# relogio e preservando todo o resto da linha.
#
# Historia, porque as duas versoes anteriores estavam erradas de formas
# diferentes e uteis:
#  - ate a rodada 16 isto era uma receita em prosa ("ignore as linhas com
#    `tempo=` e `s]`") que nao cobria a terceira forma do relogio (R16-06);
#  - a primeira versao em script declarava `#!/bin/sh` e usava substituicao de
#    processo, entao dava erro de sintaxe sob dash, que e o /bin/sh desta
#    maquina. Eu a testei com `bash script` e nunca com `./script`, que e o que
#    o README manda rodar (rodada 17, R17-06). Quarta recorrencia de "medir a
#    invocacao errada" neste projeto;
#  - a mesma versao apagava LINHAS INTEIRAS, entao engolia
#    "[ 50.5s] 300/300 raizes em 16 processos" junto com o relogio, e uma
#    contagem de raizes que mudasse passaria invisivel (R17-07).
#
# Agora o relogio e SUBSTITUIDO por <clock> no lugar, a linha sobrevive, e
# qualquer outra diferenca aparece.
#
# Uso: ./diff_ignoring_clock.sh novo.log out/comitado.log
#      ./diff_ignoring_clock.sh --selftest

norm() {
    sed -E -e 's/tempo=[0-9]+\.[0-9]+s/tempo=<clock>/g' \
           -e 's/\[[[:space:]]*[0-9]+\.[0-9]+s\]/[<clock>]/g' \
           -e 's/([0-9]+ contagens,) [0-9]+\.[0-9]+s/\1 <clock>/g' "$1"
}

if [ "$1" = "--selftest" ]; then
    ref=out/experiment_gate_richardson.log
    tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
    fail=0

    # 1. o log normalizado tem o MESMO numero de linhas que o original.
    #    (a versao anterior comparava norm(ref) com norm(ref), que e verdadeiro
    #    por construcao e nao testava nada, R18-08. Esta versao morde: e
    #    exatamente a propriedade que a versao que apagava linhas violava.)
    norm "$ref" > "$tmp/a"
    if [ "$(wc -l < "$tmp/a")" != "$(wc -l < "$ref")" ]; then
        echo "BROKEN: a normalizacao mudou a contagem de linhas"; fail=1
    fi

    # 2. o relogio E neutralizado, nas TRES formas
    sed -E -e 's/tempo=50\.5s/tempo=99.9s/' -e 's/\[  50\.5s\]/[  99.9s]/' \
           -e 's/125 contagens, 9\.7s/125 contagens, 99.9s/' "$ref" > "$tmp/clock"
    norm "$tmp/clock" > "$tmp/c"
    diff -q "$tmp/a" "$tmp/c" >/dev/null || { echo "BROKEN: o relogio nao foi neutralizado"; fail=1; }

    # 3. CONTROLE NEGATIVO: o que NAO e relogio tem de sobreviver.
    #    inclui a contagem de raizes, que fica na mesma linha do relogio e que a
    #    versao anterior apagava junto (R17-07).
    for probe in 's/media=0.62387/media=0.99999/' 's|300/300 raizes|299/300 raizes|' \
                 's/n=300/n=299/'; do
        sed -E "$probe" "$ref" > "$tmp/mut"; norm "$tmp/mut" > "$tmp/d"
        if diff -q "$tmp/a" "$tmp/d" >/dev/null; then
            echo "BROKEN: o filtro engoliu a mudanca [$probe]"; fail=1
        fi
    done

    [ $fail -eq 0 ] && echo "selftest OK: relogio neutralizado nas 3 formas, 3 controles negativos acusam"
    exit $fail
fi

[ $# -eq 2 ] || { echo "uso: $0 novo.log comitado.log" >&2; exit 2; }
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
norm "$1" > "$tmp/a"; norm "$2" > "$tmp/b"
diff "$tmp/a" "$tmp/b"
