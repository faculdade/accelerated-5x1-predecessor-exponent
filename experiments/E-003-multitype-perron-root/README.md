# E-003; the multitype Perron root of the arithmetic reverse tree (H-004)

**Read this first: the result this experiment was written to establish
did not survive review, and what replaced it is a negative result about
the method.** The first version of this README reported that the
multitype Perron root reproduces the scalar annealed exponent to twelve
decimal places, and concluded that the sibling congruence does not shift
the exponent. A `research-critic` round on 2026-08-29 showed the
agreement is an algebraic identity the construction forces, and that
every apparent deviation from it was a sampling artifact. All three
critical findings were verified independently before this rewrite.

Related hypotheses: H-004 (closed inconclusive because of this), H-002
(the identity this relies on, untouched), H-003, H-009 (opened on a
false fact from here, closed refuted the same day).
Related literature: L-041, **L-043 (Michael and Volkov 2010, whose
Corollary 1 is the published form of the criterion used here)**, L-129.

## The question it was meant to settle

E-002's theory is scalar: a fertile node is assumed to have `1/q`
expected children at each exponent, giving
`rho(alpha) = q^(alpha-1)/(2^alpha-1)` with the Kontorovich-Lagarias
value as its smaller root. That averages over residues, and H-002 proved
the residues are not averaged. If the true congruence-constrained tree
had a different multitype growth rate, the empirical reading would be
compared against the wrong number.

## Why the construction cannot answer it

Write `u' = u + q^2 t`. Then `w' = (2^a u' - 1)/q = w + q t 2^a`, and
`2^a` is invertible mod `q`, so as `t` runs over `0..q-1` the child's
residue mod `q^2` runs over all `q` lifts of `w mod q` exactly once.
Note the modulus: `w' = w + q t 2^a` leaves the child's residue mod `q`
CONSTANT, not equidistributed, and it is the spread across the `q^2`
classes that the rank-one argument needs. Instantiated at `q=5` over
the fertile `(a,u)` pairs: mod `q` gives one distinct value in every
case, mod `q^2` gives five. An earlier version of this paragraph said
`mod q` and `all of Z_q`, which is false as written; critique round 4
caught it. The conclusion below is unaffected. Averaging the
matrix rows over the `q` lifts of each parent class gives a reduced
matrix of **rank one**, whose Perron root is `q^(s-1)/(2^s - 1)`: the
scalar annealed pressure, identically, as a function of `s`, for every
`q` and at every type modulus.

Measured at `s = 0.5`, `q = 5`: the matrix has **rank 4**, equal to
`d = ord_5(2)`, on a type space of 20, and its Perron root matches the
closed-form scalar value to relative `2.1e-16`.

The row normalization `c / W[pt]` is an average over sampled roots, so
it reintroduces exactly the averaging over residues that H-002 says does
not happen. **No finite residue modulus avoids this**: a typing at `q^K`
always loses the last `q`-adic digit, and that digit decides the
grandchild's fertility.

## The three numbers that were artifacts

**Every reported deviation was governed by `n mod (d*q^2)`.**
`sample_us` walks consecutive odd integers keeping fertile residues; the
qualifying residues have period `2q^3` containing `d*q^2` of them, so
the sample is a perfect residue system mod `q^3` exactly when `n` is a
multiple of `d*q^2`. For `q = 5` that divisor is 100, and every `n` used
anywhere in the original experiment (2 000, 20 000, 50 000, 100 000) is
a multiple of 100.

| q=5, sample size | result |
|------------------|--------|
| `n = 20000` (multiple of 100) | root minus alpha_- = `-4.4e-16` |
| `n = 20037` | `+5.9e-05` |
| `n = 20013` | `-2.5e-04` |
| `n = 20050` | **no root at all** |
| `n = 100000` (multiple of 100) | `-4.4e-16` |

So `robustness.py`'s whole battery re-tested one degenerate
configuration and had no power. `main()` now prints a balanced and an
unbalanced `n` side by side, with the matrix rank, so the log exhibits
the artifact instead of hiding it.

**"The Perron root stays above 1 throughout" was false.** `matrix_root`
evaluated only the two endpoints. `rho(s)` is convex here, because
negative displacements exist (`disp = a - log2(q)`, and `a` can be 1),
so `rho(s) = 1` has two roots and equal-signed endpoints hide both.
Scanning the interior: the minimum is `0.7283` at `s = 0.50` for
`q = 11`, with 40 of 51 grid points below 1; re-bracketing recovers
roots within `5e-5` of each `alpha_-` for `q = 11, 23, 31`. Fixed: the
search now scans a grid and brackets on `[lo, argmin]`.

**The conjecture this experiment claimed to refute was never tested.**
It rested on `q = 11` and `q = 13` behaving oppositely, and both
behaviours were the artifacts above.

## What the criterion is, and whose it is

`rho(s) = 1` at the smaller root is the correct criterion, and it is
published: **Michael and Volkov 2010, Corollary 1**,

    lim_{t->inf} log E[Z(e^-t)] / t = min{ s in D : rho(s) = 1 }

The earlier attribution here, to Menshikov, Petritis and Volkov 2007,
was wrong. Their Theorem 2 gives `lambda = inf_{s>=0} rho(s)` as a
**finiteness dichotomy** for `Z(x)`, explicitly declining the critical
case `lambda = 1`, and their `m(s)` is `b x b` with exactly one child
per colour, not an unbounded sibling sum.

Corollary 1, `annealed_exact.py` and this experiment are all
**annealed**. The quenched exponent, which is what the tree realizes and
what E-001 and E-002 measure, is untouched by any of them.

## Files

| file | what it does |
|------|--------------|
| `type_closure.py` | the working version, now with the interior scan and the rank diagnostic |
| `multitype_perron.py` | a documented dead end, kept. **Carries a latent bug**: it filters residues mod `q^K` by `r % 2 == 1`, but `q^K` is odd, so a residue mod `q^K` says nothing about the parity of `u`. `build_matrix` would produce a wrong matrix if reached; it never is, because the same-modulus closure check refuses first |
| `robustness.py` | the stress tests. **They have no power**: every `n` they use is a multiple of 100 for `q = 5`, so all four re-test the same balanced configuration |

## What actually stands

**Graded closure, and it is a proof rather than a test.** If
`u = u' mod q^(K+1)` then `2^a u - 1 = 2^a u' - 1 mod q^(K+1)`, hence
`w = w' mod q^K`. One line, every `K`, every `q`. The original version
reported it as "tested at K = 1, 2, 3", and three of those cells were
vacuous: for `q = 13, 23, 31` at `K = 3` no residue class recurred in
the sample, so the check had no opportunity to fail and still printed
CLOSES.

**Fertility closure**: `u mod q^2` determines which sibling positions
are sterile. Elementary, and it is what H-002 already proves.

**Sterile children are handled correctly.** Dropping them from the
matrix is right, not a bug: an absorbing type contributes a zero row,
which cannot change a non-negative matrix's spectral radius, and the
sterile population at each level is a bounded multiple of the fertile
one.

## Reproduce

```
python3 type_closure.py 5 7 11 13 23 31   # balanced vs unbalanced n, with rank
python3 multitype_perron.py 3 5 7         # the dead end, kept for the record
python3 robustness.py                     # kept, but see the caveat above
```

## Evidence (Rule 9a)

```
Command:      python3 type_closure.py 5 7 11 13 23 31
              python3 multitype_perron.py 3 5 7
              python3 robustness.py
Commit:       d65c715 (o commit que REGISTROU estes logs). A rodada 14
              (R14-02) achou aqui `465c0bb`, que é um commit de
              aquisição de literatura: ele tocou HYPOTHESES.md e
              literature/INDEX.md e não este experimento. As
              execuções foram feitas na árvore de trabalho; o campo
              passa a nomear o commit em que a saída entrou no
              repositório, que é verificável com git log.
Date:         2026-08-29
Environment:  Linux 7.0.0-30-generic, Intel Core 7 240H (16 threads),
              62 GiB RAM; Python 3.12.3, numpy 2.5.1, scipy 1.18.0
Parallelism:  robustness.py is parallel: 18 independent units (per q,
              three sample sizes, five truncations, one large-u window),
              one per process via multiprocessing.Pool, chunksize=1 so
              the expensive (b) units do not stack on one worker.
              Measured 2026-08-29: 64.4 s serial to 20.4 s on 16
              workers, output byte-identical to the serial run (diff
              empty). sample_from and root_of carry no randomness, so a
              worker regenerating its own sample gets the same numbers
              regardless of how the work was split.
              type_closure.py e PARALELO sobre q desde a rodada 15,
              uma unidade por processo: 288.2 s para 194.9 s no comando
              declarado, saida identica (diff vazio). O ganho e modesto,
              1.48x, porque os seis valores de q tem custos muito
              diferentes e o maior domina; e uma varredura desbalanceada,
              e isso esta dito aqui em vez de apresentado como 6x.
              multitype_perron.py (0.17 s) segue serial pela razao
              medida: uma ordem de grandeza abaixo da linha de um minuto.

              TERCEIRA razao errada neste campo, e a mais instrutiva. A
              rodada 9 escreveu "type_closure.py (4.4 s) ... uma ordem de
              grandeza abaixo da linha" e a rodada 15 (R15-01) mostrou que
              a medicao foi da invocacao ERRADA: 4.4 s e o script SEM
              ARGUMENTOS, enquanto este bloco declara
              `type_closure.py 5 7 11 13 23 31`, que leva 288 s, cinco
              vezes ACIMA da linha. Medir o comando declarado, nao o
              script.
              Two wrong reasons preceded this, and both are worth
              recording. The first was "a one-off diagnostic whose
              result is now negative", which is not one of the two
              reasons the taxonomy admits, and critique round 4 refused
              it exactly as rounds 2 and 3 refused the equivalent excuse
              in E-001. The second was written while fixing the first:
              "the longest run is about 57 s, under the one-minute
              line". That number was in the block already and had never
              been measured. Measuring it gave 64.4 s, ABOVE the line,
              so the repair had replaced an inadmissible reason with a
              false one. Parallelizing was the answer both times.
Exit:         0
Output:       out/type_closure_all_q.log
Output:       out/multitype_perron.log
Output:       out/robustness.log
Nao reproduz: out/type_closure.log e out/type_closure_more_q.log sao
              saidas da versao ANTERIOR do type_closure.py, de antes da
              reescrita que a critica de fechamento da H-004 forcou. O
              formato delas ("matrix (jmax=6): ...") nao existe mais no
              script, entao NENHUM comando atual as reproduz, e isso e
              declarado aqui em vez de deixado para o revisor descobrir
              (rodada 14, R14-07). Ficam no repositorio porque registram
              o resultado que foi retratado, e a retratacao esta na
              prosa acima; nao ficam como evidencia reproduzivel.
              Os TRES comandos acima foram verificados nesta rodada:
              cada um reproduz o seu log bit a bit.
Checked:      research-critic round de FECHAMENTO DE HIPÓTESE (H-004),
              2026-08-29, REGISTRADA em notes/H-004.md (Regra 8: a
              crítica de um fechamento de hipótese mora na nota da
              hipótese), com o transcrito bruto em
              notes/codex-loop/log-B-01.txt. NÃO é uma rodada do loop do
              paper e por isso não aparece no CRITIQUE.md, que é
              específico do paper (Regra 8). A rodada 12 (R12-16)
              apontou que citar a contagem sem dizer onde ela mora
              deixava a referência solta: 3 critical, 4 major,
              5 moderate, 7 minor, overturning the headline result. All
              three critical findings were then verified independently by
              the producer before this rewrite.
              CORRECTION to the earlier evidence block, which claimed the
              q=5 root was checked against annealed_exact.py "by a
              completely different route, no matrix and no eigenvalues":
              that is false. The twelve-decimal comparison was against
              alpha_minus(), the same bisection on q^(a-1) = 2^a - 1
              duplicated verbatim in both files. There was no independent
              check, and its absence is why the tautology survived.
```

## Stale artifacts, marked rather than deleted

`out/type_closure.log` and `out/type_closure_more_q.log` are outputs of
the **pre-fix** script (bisection from `s = 0.20`, before the interior
scan) and do not reproduce from the code as committed. They are kept
because the H-009 episode is only legible with them, and they must not
be read as current results. `out/type_closure_all_q.log` is current.

`out/robustness.log` reproduz BIT A BIT, medido na rodada 15 (R15-05).
A redação anterior dizia que ele reproduzia "com diferença de último bit
do resolvedor, e não byte a byte, um artefato de threading do LAPACK",
e sugeria fixar `OMP_NUM_THREADS`. Isso não foi medido quando escrito, e
não é o caso: `diff` contra uma execução nova sai vazio.

## Sondas de verificação da rodada 3 (`critique-2026-08-29-verification/`)

Seis scripts de uso único, escritos durante uma rodada de crítica para checar
o fechamento de tipos por amostragem. **Não são parte do caminho de
reprodução do paper** e não aparecem em nenhum bloco `Command:`.

Estado, medido em 2026-08-30 ao preparar o depósito de reprodutibilidade:

- `probe3.py`, `probe4.py`, `probe5.py` rodam e saem com 0.
- `probe.py` e `probe2.py` **não rodam mais**: esperam
  `r, nty, info = matrix_root(...)`, e o `matrix_root` do `type_closure.py`
  passou a devolver um único float. Deriva de assinatura, não defeito de
  resultado; o que eles verificaram na época segue registrado no
  `CRITIQUE.md` da rodada 3.
- `probe6.py` passa de 180 s e não foi cronometrado até o fim.

Os seis carregavam um caminho absoluto para o home do autor, o que os tornava
inexecutáveis fora desta máquina; isso foi corrigido para caminho relativo em
2026-08-30. **A pasta foi deliberadamente deixada de fora do depósito
público** (`github.com/faculdade/accelerated-5x1-predecessor-exponent`): um depósito ligado a DOI
não deve embarcar código não declarado do qual um terço não roda contra a
versão atual do módulo que ele importa.
