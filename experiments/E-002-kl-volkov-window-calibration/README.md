# E-002 -- calibrating the Kontorovich-Lagarias versus Volkov gate (H-001, O8)

Related hypotheses:
[H-001](../../notes/H-001.md),
[H-002](../../notes/H-002.md).
Predecessor experiment: [`E-001`](../E-001-qx1-empirical-gate/).

## The question

Kontorovich and Lagarias (arXiv:0910.1944, Theorem 8.10) predict a
counting exponent `eta_5,BP ~ 0.650919` for the reverse tree of `5x+1`.
A competing branching model of Volkov, discussed in the same paper,
predicts `~0.678`. The two differ by `Delta = 0.027081`, and the
authors write that Volkov's data "seems insufficient to discriminate
between these two predicted exponents. It would be interesting for this
problem to be investigated further." E-001 measured `0.639` with a
bootstrap interval `[0.633, 0.645]`.

Nobody had measured what that estimator does to a process tuned to a
known annealed pressure target. Doing so turns out to answer the original
question too, because the deficit against that target is larger than
`Delta`, and the way
around that is to run the same estimator on processes tuned to each
of the two disputed exponents and see which reading the arithmetic tree
matches.

## What is here

| file | what it does |
|------|--------------|
| `tree_counts.c` | the enumerator: the arithmetic tree and the six matched stochastic controls, one code path |
| `annealed_brute.py` | the brute-force double sum behind the closed form, exact rationals, one negative-control probe per link |
| `cycle_membership.py` | which sampled roots reach a cycle, for BOTH root samples, each replica checked against its source |
| `paired_bootstrap.py` | what the band and the separation are worth; why pairing does not cancel here |
| `decade_sweep.py` | the comparison at all five checkpoint decades, with the gap's interval and the band's monotonicity probability |
| `basin_split.py` | the reading split into roots known to reach a cycle and the rest |
| `referee_tables.py` | the paper's two tables, generated from the committed count files rather than typed |
| `validate_vs_python.py` | byte-for-byte check of the C against the E-001 Python enumerator |
| `annealed_exact.py` | closed form for the exact annealed counting function of the model |
| `check_mean_vs_annealed.py` | the simulator reproduces that closed form |
| `compare_modes.py` | count distributions of the modes side by side |
| `within_root_spread.py` | separates across-root from within-root fluctuation |
| `cyc_vs_cycq.py` | checks that the integer and real-valued recursions agree |
| `buffer_squeeze.py` | measures the observed sensitivity of the truncation extrapolation to how much buffer headroom it has |
| `analyze.py` | the E-001 estimator, per-decade slopes, deficit against `alpha_-` |
| `summary.py` | the comparison table: every process, one estimator |
| `run_deep.sh` | the matched deep batch, checkpoints to `1e10`, buffers to `1e15` |

Build and reproduce:

```
gcc -O3 -march=native -fopenmp -o tree_counts tree_counts.c -lm
python3 validate_vs_python.py          # must print VALIDATION PASSED
python3 annealed_exact.py 5
./tree_counts --q 5       --roots 300 --cp 4 8 --buf 9 13 --out data/q5_arith_b13.txt
./tree_counts --q 5 --cyc --roots 300 --cp 4 8 --buf 9 13 --out data/q5_cyc_b13.txt
./tree_counts --q 5 --iid --roots 300 --cp 4 8 --buf 9 13 --out data/q5_iid_b13.txt
./tree_counts --q 5 --cycq 5.00000 --roots 300 --cp 4 8 --buf 9 13 --out data/q5_cycq500_b13.txt
./tree_counts --q 5 --cycq 5.05398 --roots 300 --cp 4 8 --buf 9 13 --out data/q5_cycq505_b13.txt
./tree_counts --q 5 --iidq 5.00000 --roots 300 --cp 4 8 --buf 9 13 --out data/q5_iidq500_b13.txt
./tree_counts --q 5 --iidq 5.05398 --roots 300 --cp 4 8 --buf 9 13 --out data/q5_iidq505_b13.txt
python3 summary.py                     # the comparison table
./run_deep.sh                          # the matched deep batch, tens of minutes
python3 summary.py b15 10              # the same table at decade 1e9->1e10
python3 buffer_squeeze.py data/q5_arith_b17.txt
```

## The five modes

They share one code path. The branch class of a node is either the true
residue or a draw, and nothing else differs:

- `arith`: `r = u mod q`, the real tree.
- `iid`: `r` drawn uniformly at every node. This is the branching random
  walk whose annealed pressure is `q^(alpha-1)/(2^alpha-1)`, so its
  annealed pressure root is `alpha_-(q) = 0.650919` at `q = 5`; whether its
  REALIZED counting exponent equals that root is assumed by the paper, not proved.
- `cyc`: the first sibling's class is drawn, and successive siblings
  advance by `c = ((2^d-1)/q) mod q`, which is what the arithmetic tree
  does exactly (H-002).
- `cycq qval`: the `cyc` structure with the value denominator replaced by
  a real `qval`, so the exponent becomes tunable. It solves
  `qval^alpha = q(2^alpha - 1)`: `qval = 5.00000` gives 0.650919 and
  `qval = 5.05398` gives 0.678.
- `iidq qval`: to `cycq` what `iid` is to `cyc`. Same tunable `qval`, but
  each child draws its own residue instead of advancing cyclically from
  one drawn at the parent. Added after an external referee observed that
  the 0.678 target had a single construction, so a spread measured on the
  0.650919 family was being carried to a target where it had never been
  checked. `cycq` and `iidq` differ in exactly one line of `dfs_root_q`.

The six controls of the paper are these five modes minus `arith`, taken
at the two targets: `cyc`, `cycq 5.00000`, `iid` and `iidq 5.00000` built
to 0.650919, and `cycq 5.05398` and `iidq 5.05398` built to 0.678.

Roots are fertile by construction in every mode. Getting this wrong
was a real error in the first pass here: the arithmetic roots are
sampled with `u mod q != 0` and are therefore always fertile, so a
control that drew the root residue from `{0..q-1}` killed one tree in
`q` outright and read `0.484` instead of `0.612`.

## The exact annealed count

For each integer `n >= 1`, a node's expected number of children at
exponent exactly `n` is `1/q`: the child exists iff `2^n r == 1 (mod q)`,
that is `r == 2^(-n)`, one residue class out of `q`. So the offspring
intensity is `(1/q) sum_{n>=1} delta_n`, a level-`k` node reached by
exponents `a_1..a_k` sits at value ratio `2^A/q^k` with `A = sum a_i`,
and

```text
E[# level-k nodes with sum a_i = A] = q^(-k) C(A-1, k-1).
```

Counting those with `2^A/q^k <= 10^t` means `A <= N_k(t)` with
`N_k(t) = floor((t + k log10 q)/log10 2)`, and the hockey-stick identity
`sum_{A=k}^{N} C(A-1,k-1) = C(N,k)` collapses the inner sum:

```text
M(t) := E[N(u0 * 10^t)] = sum_{k>=1} C(N_k(t), k) / q^k .
```

Checked by `annealed_brute.py` against the brute-force double sum for
`q = 3, 5, 7` and `t = 1..4`, in exact rational arithmetic, level by
level (`out/annealed_brute.log`), and against the simulator's mean count
by `check_mean_vs_annealed.py`. The brute-force route rebuilds the inner
count by dynamic programming over the exponent tuples, so it uses neither
the binomial coefficient nor the hockey-stick identity that produce the
closed form. Until round 16 this sentence had no implementation behind
its first half (R16-02).

This matters because it settles a question the fitted extrapolations
could not. The annealed local slope reaches `0.6517` at `t = 3` and
`0.65079` at `t = 4`, against `alpha_-(5) = 0.650919`. The annealed side
of the model has essentially no finite-window bias at the scales E-001
worked in. Whatever bias the estimator has is therefore a
quenched-versus-annealed lag, the log-slope of one realization trailing
the log-slope of the mean, and not a correction exponent that could be
fitted away.

## Result, part 1: the estimator has a bias larger than the thing it measures

Standard E-001 window, `1e5..1e8`, 300 roots, truncation extrapolated to
infinite buffer by Aitken, identical estimator in every mode:

| mode | estimator | sd of log10 N(1e8) | pressure target |
|------|-----------|--------------------|---------------|
| iid | 0.6131 | 0.8014 | 0.650919 |
| cyc | 0.6294 | 0.6657 | 0.650919 |
| arith | 0.6382 | 0.5942 | disputed |

Proveniência das duas colunas, que a rodada 14 (R14-05) apontou como
não registrada: o estimador vem de `out/summary_b13.log` e o desvio de
`out/compare_modes_b13.log`. O segundo NÃO existia; a corrida tinha sido
feita e o log nunca salvo. Salvo agora, e as três figuras conferem.

The estimator reads 0.038 below the pressure target of a process tuned to
it. Calling that a bias requires the realized-exponent assumption, which the
paper states and does not prove.
That is larger than `Delta = 0.027`. So the raw reading cannot be
compared against a theoretical prediction at all, which is what E-001
and H-001 did, and it is also why adding the bias back by hand is not
licensed: the bias itself depends on how much the process fluctuates,
and the three rows above have visibly different fluctuation.

## Result, part 2: compare readings, not a reading against a prediction

The fix is to stop comparing a window reading against a theoretical value.
Run the same estimator on a process tuned to pressure target 0.650919 and
on one tuned to pressure target 0.678, and see which reading the
arithmetic tree matches. Mode `cycq` supplies both: same branching, same
sibling congruence, same roots, same window, same buffers, with only the
value denominator changed, `qval^alpha = q(2^alpha - 1)`.

From `summary.py`:

```
       process  pressure target             window estimator              decade 1e7->1e8
  cycq 5.00000       0.650919   0.63950 [0.63357,0.64647]   0.64796 [0.64426,0.65204]
  cycq 5.05398       0.678000   0.65943 [0.65290,0.66630]   0.67079 [0.66649,0.67585]
           cyc       0.650919   0.62943 [0.62213,0.63650]   0.64437 [0.64067,0.64819]
           iid       0.650919   0.61308 [0.60233,0.62415]   0.64068 [0.63276,0.64962]
  iidq 5.00000       0.650919   0.60556 [0.59450,0.61633]   0.63202 [0.62334,0.64020]
  iidq 5.05398       0.678000   0.63264 [0.62091,0.64430]   0.65269 [0.64346,0.66146]
         arith       disputed   0.63824 [0.63183,0.64474]   0.64791 [0.64391,0.65241]
```

The arithmetic tree reads 0.64791 on the deepest decade this grid
supports. The four controls assigned pressure target 0.650919 read from
0.63202 to 0.64796, and the arithmetic reading sits inside that span.
The two controls assigned 0.678 read 0.67079 and 0.65269. Taking the
lower of the two, which is what the paper's conservative separation
does, leaves 0.00478 here, and that second control's interval
[0.64346,0.66146] overlaps the arithmetic [0.64391,0.65241]. On the
window estimator it reads 0.63264, below the arithmetic 0.63824.

This grid is a first pass, not the result. At this depth the target
deficit is still large and differs a lot between constructions. The
deep run below is the comparison the paper reports.

All seven rows use the `b13` grid, since the comparison needs the same
amount of buffer everywhere. More buffer moves the number down a little:
on the `b17` grid the arithmetic tree reads 0.6465 on that same decade
against 0.64791 here. The controls would move with it, and
`buffer_squeeze.py` puts the observed shift at 0.002, well under the 0.027
separating the two hypotheses.

So E-001's `0.639` was never evidence against Kontorovich-Lagarias. It
is, to three decimals, what a control assigned their pressure target
returns under that estimator.

The separation is conservative. Fluctuation of `log10 N(1e8)` on the
CALIBRATION grid (b15, `out/compare_modes.log`) runs 0.5941 for `arith`,
0.5976 for `cycq(5.000)`, 0.6713 for `cycq(5.05398)`, 0.6790 for `cyc`
and 0.7845 for `iid`. Until round 21 this paragraph quoted 0.594 / 0.629
/ 0.721 / 0.666 / 0.801, which are the b13 figures
(`out/compare_modes_b13.log`), under an argument about the b15
calibration. The direction survives the correction, which is why the
conclusion below is unchanged: `cycq(5.05398)` fluctuates more than
`arith` on both grids. The ordering of the middle three does not
survive, so no claim rests on it. The more variable of the two 0.678 controls also shows the larger observed
target deficit on this window; no general monotone relation between
variability and finite-depth deficit is assumed. The observed deficits are 0.019 and 0.011 on the window estimator. Reading
those as bias, and inferring from them where each realized exponent sits,
would need exactly the monotone relation just disclaimed, so no such
inference is drawn here. The bands still do not meet.

This is a measurement with calibrated controls, not a proof, and it
tests the exponent 0.678, not Volkov's model. That model is a complete
binary tree with a different encoding of the iterates, and it is not
implemented here.

Two systematics were checked rather than assumed. `cyc` and `cycq(5.0)`
are a matched integer/real-valued pair sharing the same target and
residue-label rule, and differ by 0.0035 over six
seeds, 1.5 standard errors, so no implementation systematic above about
0.004 separates the integer recursion from the real-valued one
(`cyc_vs_cycq.py`). Putting a floor at value 1 on the real-valued walk,
where the integer recursion bottoms out, changes counts but leaves the
slope identical to five decimals for five of the six seeds; the sixth
(seed 44) moves by one in the fifth decimal, 0.62852 to 0.62851. A
rodada 14 (R14-12) apontou que a redacao anterior dizia "identical to
five decimals" sem a ressalva, e o log que ela cita mostra a diferenca.

## Deep run: the same comparison where the observed target deficit is small

The deficit against the assigned target shrinks fast with depth. On the matched `b15` grid
(checkpoints `1e4..1e10`, buffers `1e9..1e15`, 300 roots, seven
datasets), the per-decade target deficit of the four controls that have
a buffer-squeeze log falls to

Derivada dos quatro `out/buffer_squeeze_*_b15.log`, coluna "all
buffers", como (pressure target menos leitura) em cada década. A
rodada 15 (R15-04) apontou que essa proveniência não estava escrita em
lugar nenhum; e ao escrevê-la apareceu que a linha `cyc` estava VELHA,
calculada antes da regeneração do dado na rodada 13, com dois valores
errados na quarta casa. Recalculada do log atual.

| process | pressure target | L=6.5 | L=7.5 | L=8.5 | L=9.5 |
|---------|---------------|-------|-------|-------|-------|
| cycq 5.00000 | 0.650919 | +0.0134 | +0.0033 | +0.0010 | +0.0011 |
| cycq 5.05398 | 0.678000 | +0.0153 | +0.0060 | +0.0005 | +0.0005 |
| cyc | 0.650919 | +0.0164 | +0.0061 | +0.0012 | -0.0003 |
| iid | 0.650919 | +0.0355 | +0.0194 | +0.0087 | +0.0034 |

At decade `1e9 -> 1e10` the largest target deficit among the six
controls is 0.0044:

```
       process  pressure target             window estimator             decade 1e9->1e10
  cycq 5.00000       0.650919   0.63263 [0.62557,0.64025]   0.64981 [0.64884,0.65075]
  cycq 5.05398       0.678000   0.65971 [0.65310,0.66628]   0.67748 [0.67651,0.67846]
           cyc       0.650919   0.63097 [0.62437,0.63778]   0.65122 [0.65014,0.65223]
           iid       0.650919   0.61250 [0.60102,0.62336]   0.64751 [0.64387,0.65044]
  iidq 5.00000       0.650919   0.61289 [0.59967,0.62748]   0.64738 [0.64440,0.64993]
  iidq 5.05398       0.678000   0.63309 [0.62198,0.64419]   0.67358 [0.67044,0.67632]
         arith       disputed   0.63809 [0.63051,0.64661]   0.64926 [0.64818,0.65027]
```

Four independent constructions assigned the 0.650919 target read
0.64738, 0.64751, 0.64981 and 0.65122 there. That band of 0.00384 is
how much the observed deficit against the assigned target still varies
with how much each one fluctuates. The arithmetic tree reads 0.64926,
inside it. The two constructions assigned 0.678 read 0.67358 and
0.67748.

The paper's conservative separation takes the LOWER of those two, inside
each bootstrap resample, against the arithmetic reading: 0.02432, with
95% interval [0.0208, 0.0273] (`paired_bootstrap.py`, `decade_sweep.py`,
`referee_tables.py`).

Until the second 0.678 calibrator (`iidq`) existed there was only one,
and this note reported 0.02822 with interval [0.0268, 0.0297] against a
band of 0.00371. Those are the superseded figures, kept here because the
ratio between them was WITHDRAWN in round 19 and the withdrawal is worth
reading: the band carries its own 95% interval of [0.00100, 0.00755], so
the ratio 0.02822/0.00371 = 7.6 actually ran [3.7, 28.2] and was never a
two-figure quantity. An earlier draft of this note said "ten", carried
into `main.tex` as "ten interval-widths" and caught there by the paper's
critique round, C-04, which flagged that no reading in the paper's own
numbers gives ten in any unit; corrected here at the source.

Pushing the arithmetic tree alone further, to checkpoints `1e12` and
buffers `1e17`, since it has no heavy tail to stall on:

| decade | slope | bootstrap | distance to 0.650919 |
|--------|-------|-----------|----------------------|
| 1e7 -> 1e8 | 0.6465 | [0.6425,0.6506] | 0.0044 |
| 1e8 -> 1e9 | 0.6487 | [0.6467,0.6506] | 0.0022 |
| 1e9 -> 1e10 | 0.6490 | [0.6479,0.6499] | 0.0020 |
| 1e10 -> 1e11 | 0.6506 | [0.6502,0.6510] | 0.0003 |
| 1e11 -> 1e12 | 0.6505 | [0.6503,0.6508] | 0.0004 |

Those bootstrap bands cover root resampling only. `buffer_squeeze.py`
measures the other term: redoing a well-buffered decade with only the
three buffers the deepest decades have available moves it by at most
0.002, and by 0.0003 to 0.0004 from decade `1e8 -> 1e9` onward. Read the
deep decades as `0.6505 +/- 0.002`, against `alpha_-(5) = 0.650919` and
`0.678`.

The window estimator itself saturates at 0.63778 by buffer `1e17`, so
E-001's Aitken value of 0.639 for the infinite-buffer limit was right,
and, under the realized-exponent assumption, the remaining deficit against
0.6509 is attributable mainly to the finite measurement window rather than
to truncation.

## Notes

- `q = 7` has `d = ord_7(2) = 3 < 6`, so only the residues in `<2> =
  {1,2,4}` are fertile: four classes out of seven are sterile, not one.
  The pressure equation is unaffected, since the expected child count per
  exponent is still `1/q`.
- The enumerator needs no visited set. The forward map is a function, so
  in the reverse tree every node has at most one parent, and a cycle
  member is reachable only from inside its own cycle; roots are sampled
  outside the three known cycles, and checked not to recur within the
  stated forward-step cutoff; unknown longer cycles are not ruled out.

## Evidence (Rule 9a)

`validate_vs_python.py`, which this README calls the gate that "must print
VALIDATION PASSED", could not be parsed by Python at all before 2026-08-29:
a migration note had been appended to it as bare text. It had therefore
never run in this repository. Repaired and run; it passes.

```
Command:      gcc -O3 -march=native -fopenmp -o tree_counts tree_counts.c -lm
              python3 validate_vs_python.py
              python3 annealed_exact.py 5
              python3 annealed_brute.py
              python3 cycle_membership.py
              python3 paired_bootstrap.py
              python3 decade_sweep.py
              python3 basin_split.py
              python3 referee_tables.py
              python3 annealed_brute.py --selftest \
                > out/annealed_brute_selftest.log
              python3 cycle_membership.py --verify-c \
                > out/cycle_membership_verify_c.log
              python3 check_mean_vs_annealed.py
              python3 summary.py
              python3 summary.py b15 10
              python3 summary.py b15 8
              python3 compare_modes.py data/q5_{arith,cyc,cycq500,cycq505,iid}_b15.txt
              python3 within_root_spread.py
              python3 cyc_vs_cycq.py
              python3 buffer_squeeze.py data/q5_arith_b17.txt
              python3 analyze.py data/q5_arith_b17.txt
              python3 analyze.py data/q5_arith_b15.txt
              for m in arith cyc cycq500 cycq505 iid; do
                python3 buffer_squeeze.py data/q5_${m}_b15.txt \
                  > out/buffer_squeeze_${m}_b15.log; done
              python3 compare_modes.py data/q5_arith_b13.txt \
                data/q5_cyc_b13.txt data/q5_cycq500_b13.txt \
                data/q5_cycq505_b13.txt data/q5_iid_b13.txt \
                > out/compare_modes_b13.log
Commit:       ba89add (o commit que REGISTROU estes logs). A rodada 14
              (R14-02) achou aqui `465c0bb`, que é um commit de
              aquisição de literatura: ele tocou HYPOTHESES.md e
              literature/INDEX.md e não este experimento. As
              execuções foram feitas na árvore de trabalho; o campo
              passa a nomear o commit em que a saída entrou no
              repositório, que é verificável com git log.
Date:         2026-08-29
Environment:  Linux 7.0.0-30-generic, Intel Core 7 240H (16 threads), 62 GiB RAM; Python 3.12.3, numpy 2.5.1, scipy 1.18.0; gcc with -O3 -march=native -fopenmp
Parallelism:  the C enumerator is OpenMP over roots and was built with
              -fopenmp on 16 threads. The Python scripts split in two
              groups, and the reason differs by group (Rule 9c).
              summary.py and buffer_squeeze.py only post-process a
              committed count file, and annealed_exact.py evaluates a
              closed form with no input at all; all three finish in
              seconds, so they are serial for the obvious reason.
              within_root_spread.py, cyc_vs_cycq.py, validate_vs_python.py
              and check_mean_vs_annealed.py do NOT read a committed count
              file: each one drives tree_counts through subprocess.run,
              and takes 17 to 100 s. They are serial DELIBERATELY, because
              every tree_counts call is itself OpenMP over all 16 threads;
              issuing those calls from a process pool would oversubscribe
              the machine and run slower than the serial driver.
              annealed_brute.py, paired_bootstrap.py, decade_sweep.py and
              basin_split.py are serial because each finishes in a few
              seconds over committed count files. cycle_membership.py is
              serial at about 9 s, which includes its audit pass at twice
              the step cutoff (R19-01).
              An earlier version of this field claimed all the Python
              scripts post-process a committed file, which is false for
              those four (round 16, R16-03).
Exit:         0 for every command above
Output:       out/validate_vs_python.log
Output:       out/annealed_exact_q5.log
Output:       out/annealed_brute.log
Output:       out/cycle_membership.log
Output:       out/paired_bootstrap.log
Output:       out/decade_sweep.log
Output:       out/basin_split.log
Output:       out/referee_tables.log
Output:       out/annealed_brute_selftest.log (controle negativo, uma sonda
              por elo da cadeia)
Output:       out/cycle_membership_verify_c.log (confere a replica do
              splitmix64 contra um probe compilado, ALEM da checagem contra a
              coluna de raizes dos arquivos comitados, que roda sempre)
Output:       out/check_mean_vs_annealed.log
Output:       out/summary_b13.log
Output:       out/summary_b15_d10.log
Output:       out/summary_b15_d8.log
Output:       out/compare_modes.log
Output:       out/compare_modes_b13.log
Output:       out/within_root_spread.log
Output:       out/cyc_vs_cycq.log
Output:       out/buffer_squeeze_arith_b17.log
Output:       out/analyze_arith_b17.log
Output:       out/analyze_arith_b15.log
Output:       out/buffer_squeeze_arith_b15.log
Output:       out/buffer_squeeze_cyc_b15.log
Output:       out/buffer_squeeze_cycq500_b15.log
Output:       out/buffer_squeeze_cycq505_b15.log
Output:       out/buffer_squeeze_iid_b15.log
Gap found:    the paper's shallow-decade comparison quoted 0.64622 and
              0.67200 at five decimals, and no committed output held
              those numbers. The four-decimal buffer-squeeze logs round
              to them, and the only five-decimal log for that decade was
              the OTHER grid (b13, which reads 0.64791 and 0.67079), so
              the quoted pair was not reproducible from what this folder
              contained. `summary.py b15 8` had been run and its log
              never saved. Run and saved on 2026-08-29: it returns
              0.64622 and 0.67200 exactly, so the paper was right and
              the evidence was missing, not the number wrong.
Checked:      producer, against the calibrated numbers in the paper's
              abstract and \S"A calibrated comparison". All reproduce:
              arithmetic tree 0.64926; the three constructions built to
              0.650919 read 0.64981, 0.65122, 0.64751, so the largest
              distance is 0.00341 and the band is 0.00371 wide; the
              construction built to 0.678 reads
              0.67748, which is 0.02822 away (paper: 0.02821; the 0.02822
              here is the difference of the two ROUNDED readings), with paired
              95% interval [0.0268, 0.0297] (paper: the same). The band-width
              ratio this block used to certify was withdrawn from the paper in
              round 19; certifying agreement with a sentence main.tex no longer
              contains is what R20-01 caught, and the external round X11 found
              the same thing had happened to the two figures just above, whose
              "(paper: ...)" tags are gone for that reason: since the second
              0.678 calibrator entered, the 0.650919 family has FOUR members,
              its band is 0.00384, and the paper's headline is the conservative
              0.02432 with [0.0208, 0.0273]. `referee_tables.py` and
              `out/referee_tables.log` are the current record. What survives
              here unchanged is the gap against the retuned relaxation ALONE,
              which `main.tex` still carries. The unconstrained
              construction reads 0.61308 (paper: 0.6131), which is below its
              pressure target by
              0.03784 (paper: 0.038). The annealed closed form converges to
              0.650919, the Kontorovich-Lagarias value itself.
```

## Regeneração do q5_cyc_b15.txt (rodada 13, R13-01)

**O arquivo comitado não reproduzia mais a partir do código.** Este é o
defeito mais caro que um depósito com DOI pode carregar, porque o DOI é
imutável.

Causa: a guarda `if (u == 0) { zero_nodes++; continue; }` em
`tree_counts.c` pula o sorteio do RNG do `MODE_CYC` oito linhas abaixo,
deslocando o fluxo de resíduos. Antes dela, `u == 0` caía adiante, `w`
estufava por underflow para 3.69e18 e o nó virava folha **por
acidente**, mas o sorteio acontecia.

O comentário da guarda afirmava que "os dados comitados não são
afetados (verificado re-rodando e diffando)". Era falso e a verificação
nunca tinha sido feita.

Provado por teste negativo executado: uma compilação com aquela linha
removida, e só ela, reproduz o arquivo antigo exatamente.

### O que foi medido, não suposto

Os ONZE arquivos `q5_*.txt` foram regenerados. Dez reproduzem bit a bit;
só o `cyc` divergia, em duas das 300 raízes (557 e 5621) e na linha
`#resid`. O `q5_arith_b17.txt`, a corrida profunda de 824 s (cronometrado na rodada 15; o 1026 s anterior foi medido com a maquina carregada), reproduz.

O dado foi regenerado a partir do código atual e as análises que
dependem dele re-rodadas:

```
                          comitado          regenerado
  cyc, janela             0.63097           0.63097     igual
  cyc, 1e9->1e10          0.65122           0.65122     igual
  cyc, 1e7->1e8           0.64482           0.64482     igual
  IC do bootstrap    [0.62435,0.63777] [0.62437,0.63778]  5a casa
  buffer-squeeze 1e5->1e6 0.6136            0.6137       4a casa
  compare_modes mean N    167.0             227.4        muda
```

Nenhum número que o paper imprimia então se move: 0.65122, 0.00141,
0.0037, 0.0282 e 0.0263 seguem idênticos, conferidos um a um contra o
`main.tex` NA DATA DESTA ENTRADA. Os que mudam não aparecem no paper,
também conferido.

Dos cinco, só 0.65122 e 0.0282 continuam no `main.tex` de hoje. 0.00141,
0.0037 e 0.0263 saíram quando o segundo calibrador de 0.678 entrou e a
família de 0.650919 passou de três para quatro construções. A lista já
tinha perdido um sexto número, 7.6, removido na rodada 20 porque saiu do
paper na 19; a rodada externa X11 achou que os outros três precisavam do
mesmo tratamento, e o padrão é sempre este: um registro que afirma
concordância com o `main.tex` envelhece toda vez que o `main.tex` muda,
e nada avisa.

Os quatro logs derivados foram regravados: `summary_b15_d10.log`,
`summary_b15_d8.log`, `buffer_squeeze_cyc_b15.log`, `compare_modes.log`.

`validate_vs_python.py` volta a passar (`VALIDATION PASSED`) com o
enumerador recompilado.
