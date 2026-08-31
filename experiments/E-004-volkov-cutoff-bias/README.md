# E-004; the direction and size of the cutoff bias in Volkov's estimates

Related hypothesis: H-013 (closed-confirmed). Related literature: L-002
(Volkov 2006, the primary source, read on rendered pages because its
text layer is corrupt).

## The question

Volkov's conjecture that the non-escaping count grows like `K^0.678` is
supported in his paper by four numerical estimates, computed with an
artificial threshold `T = 4.2e8` standing for infinity:

```text
beta_hat(1e3)=0.741  beta_hat(1e4)=0.715
beta_hat(1e5)=0.694  beta_hat(1e6)=0.680
```

He stops there and writes that larger `m` gets too close to `T`.
Kontorovich and Lagarias cite this data as insufficient to discriminate
between `0.650919` and `0.678`, which is the opening of the paper this
experiment supports. Nobody had asked in which DIRECTION the threshold
moves the estimate.

## The answer, and it is one-sided

A seed that exceeds `T` and would later return to a cycle is counted as
escaping. No genuinely escaping seed is ever counted as returning.
Therefore `q(m,T) <= Qtilde(m)` for every `T`, and the bias in
`log q / log m` is **downward**. Raising the threshold can only raise
the estimate.

Measured, at `m = 1e6`:

| threshold T | "good" count | beta_hat |
|-------------|--------------|----------|
| 4.2e8 (his) | 12250 | 0.68136 |
| 1e9  | 12621 | 0.68352 |
| 1e10 | 13051 | 0.68594 |
| 1e12 | 13232 | 0.68694 |
| 1e16 | 13322 | 0.68743 |

The correction is `+0.0061`, and it moves the estimate **away** from
`0.650919` and **above** Volkov's own `0.678`.

## What this does NOT establish

The finite-`m` correction is a separate matter and its sign is not
determined. Writing `Qtilde(m) = m^beta L(m)`, the estimator's error is
`log L(m)/log m`, and nothing in the conjecture fixes the sign of
`log L`. Volkov's own sequence falls, which suggests a downward
pre-asymptotic drift, but the threshold error rides on it too, and not
evenly: raising T to 1e16 moves the four values by +0.00000, +0.00015,
+0.00212 and +0.00607 as m runs from 1e3 to 1e6, so the error grows with
m and steepens the fall instead of shifting it. The sequence therefore
does not isolate the finite-`m` effect. An earlier version of this
paragraph said the sequence carries the threshold error at every `m`,
which critique round 5 instantiated and refuted: at m = 1e3 the error is
exactly zero.

The honest conclusion is narrower than "the data points the other way":
his four numbers are consistent with `0.678` and cannot discriminate it
from `0.650919`, and this experiment makes that insufficiency
quantitative rather than asserted.

## The map, and a version of this file that had it wrong

Volkov's equation (1.1), page 2, read from a rendered image:

```text
M(x) = 5x + 1   if x is odd
M(x) = x / 2    if x is even
```

**The first version of this experiment implemented `(5x+1)/2` instead**,
which is Kontorovich and Lagarias's `T_5`, a different map. Volkov's
orbit passes through the peak `5x+1`; that version recorded half of it,
so the threshold test never saw the peak, and every seed whose peak
landed in `(T, 2T]` was classified wrongly. Found by critique round 3.

The error is visible in the reproduction quality. With the wrong map the
third estimate came out `0.69499`, rounding to `0.695` against his
printed `0.694`, and the fourth missed by `0.0035`. With the map as
printed the third is `0.69423`, rounding to `0.694` exactly, and the
fourth misses by `0.0014`. The residual gap at `m = 1e6` means some
detail of his implementation is still not recoverable from the text.

The cycles are computed here rather than copied, because Volkov's own
text prints the first one as `1 -> 6 -> 3 -> 1` while `M(3) = 16`; the
printed cycle is incomplete. Computed, it is
`1 -> 6 -> 3 -> 16 -> 8 -> 4 -> 2 -> 1`, and the union of the three is
27 values.

**Why the wrong version looked validated.** At `T = 1e16` both
implementations return `13322`, because at that depth the factor-of-two
difference no longer changes any classification. That agreement was
taken as confirmation against an independent count. It confirmed only
the endpoint the two maps share.

## Reproduce

```
gcc -O2 -o volkov_count volkov_count.c -lm
./volkov_count 1000000 420000000
./volkov_count 1000000 10000000000000000
```

## Evidence (Rule 9a)

```
Command:      ./run.sh > out/volkov_cutoff.log
              (reproduz o log INTEIRO, cabecalho incluido, bit a bit)
Commit:       c2aac41 (o commit que REGISTROU estes logs). A rodada 14
              (R14-06) achou aqui `e6bd652`, que e o commit da rodada 2 de
              critica e nao a execucao
Date:         2026-08-29
Environment:  Linux 7.0.0-30-generic, Intel Core 7 240H (16 threads),
              62 GiB RAM; gcc with -O2
Parallelism:  serial, pela razao MEDIDA: o `./run.sh` declarado acima,
              que faz DEZESSETE invocacoes do programa, leva 3.57 s
              (cronometrado na rodada 15), contra a linha de "roughly a
              minute of projected wall
              clock" que .claude/rules/experiments.md fixa para
              paralelizar primeiro (Rule 9c).

              A rodada 15 (R15-06) apontou que o numero anterior, 1.87 s,
              media OITO execucoes enquanto o Command declarado faz
              dezessete. Mesmo erro de forma que o R15-01 no E-003: medir
              coisa diferente da que o campo declara. A conclusao nao
              muda, porque 3.57 s segue muito abaixo da linha, mas o
              numero agora e o do comando. A forma e uma particao de
              faixa sobre sementes e paralelizaria trivialmente; o que a
              mantem serial e o tempo, e agora ele esta medido. A
              redacao anterior dizia "under a minute" sem cronometrar
              nada, e a rodada 11 (R11-07) apontou isso.
Exit:         0
Output:       out/volkov_cutoff.log
Checked:      critique round 3, which found the map was wrong and
              reimplemented Volkov's rule independently; its numbers and
              this file's now agree. Also checked by the producer against
              the printed values on page 19 of L-002 and the map on page
              2, both read from rendered images because that PDF's text
              layer is corrupt. The cycles are computed by the program
              itself from the seeds 1, 13 and 17, not transcribed.
```

## Fatos lidos na fonte primária (L-002, preprint do autor)

Lidos em imagens renderizadas a 130 dpi, porque o PDF usa fontes Type 3
do Ghostscript 6.52 e a camada de texto sai como lixo. Página do PDF,
não do periódico.

- p. 2, eq. (1.1): `M(x) = 5x+1` se `x` é ímpar, `x/2` se `x` é par.
  Não é o mapa acelerado. A rodada 3 achou o crítico R3-07 justamente
  porque a versão anterior deste programa usava `(5x+1)/2`.
- p. 2: os três ciclos impressos, com o primeiro escrito como
  `1 → 6 → 3 → 1`. Está incompleto: `M(3) = 16`, então o ciclo é
  `1 → 6 → 3 → 16 → 8 → 4 → 2 → 1`. Por isso o programa calcula os
  ciclos a partir das sementes 1, 13 e 17 em vez de transcrevê-los.
- p. 2: "we conjecture that there are only three cycles mentioned
  above". A finitude é conjectura dele, não teorema.
- p. 19, Conjectura 1: `lim (log Q̃(K) / log K) = β_M`. O limite
  determinístico é igualado ao valor do modelo, o que é a razão de o
  paper poder ler os dois valores contra a mesma quantidade.
- p. 19, seção 5.1: limiar `T = 4.2 × 10^8`; `k` é *good* se
  `M^(n)(k) ∈ {1, 13, 17}` para algum `n`; `β̂(m) = log q(m) / log m`.
- p. 19: os quatro valores impressos são `0.741`, `0.715`, `0.694` e
  `0.680`. Este programa devolve os três primeiros e `0.681` no quarto.
  A diferença de uma unidade na última casa em `m = 10^6` está
  declarada no paper, e não foi explicada: o texto não diz em que
  precisão aritmética ele contou.

## run.sh, e por que ele passou a existir (rodada 14, R14-06)

O bloco de evidência trazia como comando
`./volkov_count <m> <T> for m in 1e3..1e6 and T in 4.2e8..1e16`. Isso é
um **esquema**, não um comando: não roda como está escrito. E o
`out/volkov_cutoff.log` carregava linhas de cabeçalho que **eu escrevi à
mão** e que o programa não emite. Um revisor que clonasse não teria como
reproduzir o arquivo.

`run.sh` produz o log inteiro por máquina, cabeçalho incluído, e foi
verificado reproduzindo bit a bit o arquivo comitado.

Uma armadilha encontrada ao escrevê-lo, que vale para todo o depósito:
sem `LC_ALL=C` o `awk` usa vírgula decimal no locale pt_BR, e o log sai
diferente conforme a máquina de quem roda. O script fixa o locale.
