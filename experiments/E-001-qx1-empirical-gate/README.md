# E-001; Empirical statistical gate for the qx+1 generalization (H-001)

Related hypothesis: [H-001](../../notes/H-001.md)

## What was done

Measures the empirical counting exponent of the real 5x+1 reverse tree
to decide between two competing theoretical predictions from the
literature: Kontorovich-Lagarias (2009), 0.650919 (identical to our
second root of the qx+1 pressure equation, H-109), vs. Volkov, 0.678.

## Files

- `pressure_qx1.py`, `empirical_qx1_tree.py`, the Fable's original
  scripts (multitype pressure-equation verification; the first ad hoc
  empirical confirmation, cited in H-109 but never persisted in the
  repository before).
- `pilot_gate_5x1.py`, `pilot2_gate_5x1.py`; noise/bias/cost
  calibration pilots (q=5 admissibility rule, sensitivity to
  truncation bias).
- `experiment_gate_production.py`; production run, n=300, fixed
  window 1e5-1e8, buffer=5 decades (first version; documented in H-001
  as the step that revealed the "CI excludes both candidates" failure
  mode predicted by the Fable).
- `experiment_gate_richardson.py`; final version: DFS with path-max
  tracking (a single pass gives the counts at every buffer
  simultaneously, validated byte-for-byte against the old method) +
  Aitken Δ² extrapolation on the mean curve across roots + bootstrap.

## Result (original, 2026-07-17; superseded reading, kept for the record)

Aitken Δ² (buffer→∞): **0.639, 95% CI=[0.633, 0.645]**, read at the
time as excluding Volkov (0.678) with wide margin (~10+ standard
errors), with the residual gap to Kontorovich-Lagarias (0.650919)
attributed to a fixed-window pre-asymptotic (the per-decade slope panel
was still rising in the last tested decade), not to uncorrected
truncation bias.

**This exclusion claim does not hold.** H-137 (2026-08-07) and E-002
(2026-08-09) show the same estimator reading 0.038 below the assigned
pressure target of an unconstrained control, more than the 0.027081
that separates the two disputed values, so a raw 0.639 reading was
never evidence against either prediction. Calling that difference a
bias requires equality between a control's realized exponent and its
assigned target, which the paper states and does not prove; read as a
distance, it is close to what the estimator returns on a control
assigned 0.650919.

E-002 calibrates against matched synthetic controls that share the
arithmetic tree's branching and sibling-spacing law. The current
experiment uses four controls assigned target 0.650919 and two assigned
0.678. At the deepest matched checkpoint decade the arithmetic tree
reads 0.64926, inside the span [0.64738, 0.65122] of the first four,
and 0.02432 below the lower of the two assigned 0.678, with 95% joint
bootstrap interval [0.0208, 0.0273] on that difference, the smaller of
the two 0.678 readings being taken inside each resample.

Until the second 0.678 calibrator existed this paragraph reported a
band of three constructions and a difference of 0.0282 with interval
[0.0268, 0.0297], and a ratio in band-widths on top of that. The ratio
was withdrawn in round 19, when the band's own interval was measured,
and round 20 (R20-01) found it still alive in five files outside
main.tex. The 0.0282 figure survives in the paper only as the gap
against the retuned relaxation alone, never as the headline. See H-001 for
the full history, `../E-002-kl-volkov-window-calibration/` for the
calibrated experiment, and H-001 also for a necessary correction to an
earlier H-109 claim (the "1.547 vs 1.5363" Hill estimator cited there
is not confirmatory, the real standard error is ~0.45).

## Reproduce

```
python3 experiment_gate_richardson.py
```


Migration note: references `H-109`, `H-137`, out of this repo's scope, left unrenumbered.

## Evidence (Rule 9a)

Before 2026-08-29 this folder held scripts and no output at all, while the
paper quoted `0.639` from it. Two of its scripts did not even parse: a
migration note had been appended as bare text, not as a comment. Both were
repaired and the production script was run.

```
Command:      python3 experiment_gate_richardson.py
              python3 pressure_qx1.py > out/pressure_qx1.log
              python3 empirical_qx1_tree.py > out/empirical_qx1_tree.log
              python3 validate_singlepass_vs_perbuffer.py \
                > out/validate_singlepass_vs_perbuffer.log
              ./diff_ignoring_clock.sh --selftest \
                > out/diff_ignoring_clock_selftest.log
Commit:       c2aac41 (o commit que REGISTROU estes logs). A rodada 14
              (R14-02) achou aqui `465c0bb`, que é um commit de
              aquisição de literatura: ele tocou HYPOTHESES.md e
              literature/INDEX.md e não este experimento. As
              execuções foram feitas na árvore de trabalho; o campo
              passa a nomear o commit em que a saída entrou no
              repositório, que é verificável com git log.
Date:         2026-08-29
Environment:  Linux 7.0.0-30-generic, Intel Core 7 240H (16 threads), 62 GiB RAM; Python 3.12.3, numpy 2.5.1, scipy 1.18.0; gcc with -O3 -march=native -fopenmp
Parallelism:  parallel over roots, 16 worker processes via
              multiprocessing.Pool. The shape is independent trials, one
              per root, which this project classifies as
              always parallel. Wall clock fell from 311.3 s serial to
              51.0 s. Pool.map preserves input order, so the bootstrap
              resamples the same list it saw serially and every number
              is unchanged: the per-buffer curve, the Aitken value
              0.63896 and the interval [0.63283, 0.64522] all reproduce
              exactly. An earlier version of this block justified the
              serial run instead of fixing it, which critique rounds 2
              and 3 both refused (Rule 9c).
Exit:         0
Output:       out/experiment_gate_richardson.log
Output:       out/pressure_qx1.log
Output:       out/empirical_qx1_tree.log
Output:       out/validate_singlepass_vs_perbuffer.log
Output:       out/diff_ignoring_clock_selftest.log (controle negativo; a
              rodada 19 (R19-06) achou que nenhum controle negativo do
              deposito tinha execucao registrada sob a Regra 9a, entao um
              revisor clonando o repositorio nao tinha como saber que
              existiam)
Nao bit a bit: `experiment_gate_richardson.log` e
              `validate_singlepass_vs_perbuffer.log` embutem tempo de
              RELOGIO na saida, entao nao reproduzem byte a byte por
              construcao: so as linhas de tempo mudam. Declarado aqui na
              rodada 15 (R15-08) em vez de deixado para o revisor
              descobrir que um `diff` nao sai vazio. Todo NUMERO neles
              reproduz; um revisor que queira comparar roda
              `./diff_ignoring_clock.sh novo.log out/comitado.log`, que
              sai vazio para os dois. A receita em prosa que estava aqui
              ("ignore as linhas com `tempo=` e `s]`") nao cobria a
              terceira forma do relogio, `... contagens, 9.7s`, e quem a
              seguisse ao pe da letra ainda via diff (rodada 16, R16-06).
              O script tem controle negativo: `--selftest` roda tres
              mutacoes (a media por buffer, a contagem `300/300 raizes` e
              o `n=300`) e confirma que o filtro NAO engole nenhuma, alem
              de checar que a normalizacao preserva a contagem de linhas.
              Esta ultima existe porque a versao anterior do filtro
              apagava LINHAS INTEIRAS e engolia a contagem de raizes
              junto com o relogio (rodada 17, R17-07); o README dizia
              "uma media" ate a rodada 18 achar que as tres nunca
              chegaram aqui (R18-05).
              A rodada 14 (R14-07) fez esta mesma declaracao no E-003 e
              ela NAO viajou para ca, que e a Regra 8d2 outra vez
              (R15-07).
Also run:     `pressure_qx1.py` (exit 0, out/pressure_qx1.log) and
              `empirical_qx1_tree.py` (exit 0, 5.0 s, seeded
              random.seed(123), out/empirical_qx1_tree.log). Both were
              listed in this folder with no output of their own until
              2026-08-29. The first returns the smaller pressure root
              0.650918639898 for q=5 and the tail exponent 1.536290311423,
              independently of E-002, and confirms the multitype
              Perron root equals the scalar closed form to 1e-15 across
              q in {3,5,7} and k up to 4, which is the same collapse
              E-003 later found the hard way. The second supplies the
              evidence for the Hill retraction stated above: the tail
              index reads 1.547, 1.475 and 1.766 at the top 2%, 5% and
              10% of 600 roots, so the 1.547 that H-109 quoted against a
              prediction of 1.536 was the one fraction of three that
              agreed, and the spread across fractions is what the
              retraction rests on.
Bug fixed:    the buffer bucketing used
              `bisect_right(buffers, log10(pmax) + 1e-9)`. The epsilon
              pushed a path-max just under a power of ten into the next
              buffer: for `pmax = 999999999` the node was charged to
              `1e10` instead of `1e9`. Found in critique round 5
              (R5-15). Replaced by integer thresholds,
              `bisect_left(bufthr, pmax)`, which is what the C
              enumerator has always done (tree_counts.c:215). Re-run
              after the change: every published number is unchanged
              (per-buffer curve, per-decade slopes, Aitken 0.63896,
              interval [0.63283, 0.64522]); the only diff against the
              previous log is the wall-clock line. Verified by re-running,
              not by arguing the window was too narrow to matter.
Cross-check:  `validate_singlepass_vs_perbuffer.py` (exit 0, 9.7 s,
              out/validate_singlepass_vs_perbuffer.log). The paper's
              Method says the single depth-first pass with path-max
              tracking "agrees, count for count, with an independent
              buffer-by-buffer implementation on five of the 300 roots".
              Critique round 4 found that sentence had no committed run
              behind it, and that this block wrongly said nothing in the
              paper rested on `experiment_gate_production.py`. It does.
              The check now exists and compares the two ALGORITHMS,
              which is a different thing from E-002's
              out/validate_vs_python.log (the same single-pass routine
              in C and in Python): decade_counts_2d_v2, one pass, five
              buffers at once, against decade_counts, which re-runs the
              whole enumeration with a different hard bound per buffer.
              5 roots x 5 buffers x 5 checkpoints = 125 counts, all
              equal. The five roots are 2051, 8333, 1781, 3759, 9947,
              the same five the production panel uses, because both
              scripts sample with random.Random(2026) and an identical
              sample_roots.
              Written on 2026-08-29, and it failed on its first run:
              decade_counts_2d_v2 returns raw[checkpoint][buffer] and
              the comparison had the two indices transposed. The bug
              was in the new check, not in either enumerator.
Not re-run:   `experiment_gate_production.py` is now exercised through
              the cross-check above, which imports its `decade_counts`.
              Its own `main()` was not run: it computes the SUPERSEDED
              reading this README retracts above. `pilot_gate_5x1.py`
              and `pilot2_gate_5x1.py` were repaired (they did not
              parse) and now do parse, but were NOT executed, and
              nothing in the paper rests on their output. Stated rather
              than left silent, which is what Rule 9a exists to prevent.
Checked:      producer, against every number `papers/04-.../main.tex`
              quotes from this experiment. All reproduce exactly: the
              per-buffer curve 0.60049, 0.62387, 0.63261, 0.63650, 0.63801;
              the per-decade slopes 0.6021, 0.6296, 0.6432, 0.6460; the
              Aitken value 0.63896 (paper: 0.639); the bootstrap interval
              [0.63283, 0.64522] (paper: [0.633, 0.645]).
```

## Precisão da curva por buffer (rodada 9, R9-08)

O log imprime as médias por buffer com cinco casas. O paper enuncia a
razão de encolhimento dos incrementos como faixa, e no terceiro dígito
essa faixa depende de casas que o log não tem.

Das médias impressas, os incrementos são 0.02338, 0.00874, 0.00389 e
0.00151, e as razões 0.37382, 0.44508 e 0.38817. A linha
`incrementos:` do próprio log, que diferencia antes de arredondar, dá
0.00388 no terceiro e portanto razões 0.44394 e 0.38918. As duas rotas
já discordam na quarta casa, que é exatamente o que este parágrafo
afirma. Propagando o
arredondamento de mais ou menos 5e-6 em cada média, as razões ficam em
[0.3732, 0.3744], [0.4434, 0.4467] e [0.3846, 0.3918].

Ou seja: o log NÃO decide entre 0.44 e 0.45 para a razão do meio. O
paper dizia "0.37 to 0.44", que a envoltória não sustenta, e passou a
dizer "0.37 to 0.45", que ela sustenta. Uma versão futura que precise do
terceiro dígito tem de imprimir as médias com mais casas antes de
enunciar a faixa, não depois.
