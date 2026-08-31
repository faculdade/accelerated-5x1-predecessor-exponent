#!/usr/bin/env python3
"""
E-002, part 5: do the sampled roots reach a cycle, or do they escape?

The paper's bridge from the measured exponent to Volkov's quantity leans on
Conjecture 7.2, and how much weight the conjecture carries depends on how many
sampled roots sit in a known cycle's basin.

THIS PROJECT HAS TWO ROOT SAMPLES OF SIZE 300 AND THEY ARE NOT THE SAME SET.
They share 31 roots. Reporting one where the other is meant is exactly the error
critique round 17 caught (R17-01), after round 16 introduced it by writing this
script against the wrong one:

  E-001 sample : Python `random.Random(2026)` driving `sample_roots` in
                 E-001/experiment_gate_richardson.py. Backs sections 4 and 5.
  E-002 sample : the C enumerator's own splitmix64, seeded 2026, in
                 tree_counts.c `main`. Backs section 6, the calibrated
                 comparison, and every `data/q5_*` count file.

So this script reports BOTH, labelled, and refuses to run if either sampler
stops matching its source. Both checks are UNCONDITIONAL, on the declared
command, because round 18 (R18-01) found the previous version's check gated
behind a `--verify-c` flag the README never passes: re-injecting round 16's
exact mistake there produced 19 and 425-549 under the E-002 label, exit 0, with
nothing caught.

What each replica is checked against matters as much as that it is checked:

  E-001 : against the real `sample_roots`, imported from E-001 and called.
  E-002 : against the FIRST COLUMN of every committed `data/q5_*.txt`, which is
          the root list the C enumerator actually wrote when it produced the
          numbers the paper quotes.

The E-002 check deliberately does NOT compare against a re-transcription of
tree_counts.c's sampler into this file. A transcription agreeing with itself
proves nothing (CLAUDE.md Rule 8g, item 1: observe the artifact, never the
process). `--verify-c` still compiles such a probe, and it is an extra, not the
check.

Forward accelerated 5x+1: from odd u, u -> (5u+1)/2^v where 2^v || 5u+1,
iterated up to 5000 steps. Every root that reaches a cycle does so within 42
steps in both samples. The script does not take that on trust, and does not
audit it with a probe below the cutoff either, which is what R18-03 and R19-01
each got wrong: it re-runs at TWICE the cutoff and requires the same count, so
the check fails exactly when a root enters a cycle past 5000.

Serial (Rule 9c): 600 orbits at the 5000-step cutoff plus 600 more at the
10000-step audit, on integers reaching ~1100 digits, take about 9 s together,
under the one-minute line at which .claude/rules/experiments.md asks for
parallel decomposition.

Run: python3 cycle_membership.py
     python3 cycle_membership.py --verify-c   # also cross-check vs the C sampler
"""
import glob
import os
import random
import subprocess
import sys
import tempfile

CYCLE5 = {1, 3, 13, 33, 83, 17, 27, 43}
N_ROOTS, STEPS, LO, HI, Q, SEED = 300, 5000, 101, 9999, 5, 2026
M64 = (1 << 64) - 1


def sample_e001(n=N_ROOTS):
    """Python random.Random(2026); verbatim rule from
    E-001/experiment_gate_richardson.py:92. Backs sections 4 and 5."""
    rng = random.Random(SEED)
    roots, seen = [], set()
    while len(roots) < n:
        v = rng.randrange(LO, HI, 2)
        if v % Q and v not in CYCLE5 and v not in seen:
            seen.add(v)
            roots.append(v)
    return roots


def _sm64(s):
    """splitmix64, transcribed from tree_counts.c:81"""
    s = (s + 0x9E3779B97F4A7C15) & M64
    z = s
    z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & M64
    z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & M64
    return s, z ^ (z >> 31)


def sample_e002(n=N_ROOTS):
    """the C enumerator's own sampler, transcribed from tree_counts.c:365-379.
    Backs section 6 and every data/q5_* count file."""
    s, roots, seen = SEED, [], set()
    while len(roots) < n:
        s, r = _sm64(s)
        v = LO + (r % (HI - LO + 1))
        if v % 2 == 0 or v % Q == 0 or v in CYCLE5 or v in seen:
            continue
        seen.add(v)
        roots.append(v)
    return roots


C_PROBE = r"""
#include <stdio.h>
#include <stdint.h>
static inline uint64_t sm64(uint64_t *s){uint64_t z=(*s+=0x9E3779B97F4A7C15ULL);
 z=(z^(z>>30))*0xBF58476D1CE4E5B9ULL;z=(z^(z>>27))*0x94D049BB133111EBULL;return z^(z>>31);}
static int in_cycle(uint64_t v){static const uint64_t c[]={1,3,13,33,83,17,27,43};
 for(int i=0;i<8;i++) if(v==c[i]) return 1; return 0;}
int main(void){uint64_t r[300],s=2026,lo=101,hi=9999;int q=5,h=0;
 while(h<300){uint64_t v=lo+(sm64(&s)%(hi-lo+1));
  if(!(v&1))continue; if(v%q==0)continue; if(in_cycle(v))continue;
  int d=0; for(int i=0;i<h;i++) if(r[i]==v){d=1;break;} if(d)continue; r[h++]=v;}
 for(int i=0;i<300;i++) printf("%llu\n",(unsigned long long)r[i]); return 0;}
"""


def verify_against_data(quiet=False):
    """the E-002 replica against the artifact: the root column of the committed files.

    These files are what the C enumerator wrote when it produced the counts the
    paper quotes, so agreeing with them is agreeing with the real sampler."""
    here = os.path.dirname(os.path.abspath(__file__))
    files = sorted(glob.glob(os.path.join(here, "data", "q5_*.txt")))
    if not files:
        raise SystemExit("ERROR: no data/q5_*.txt to check the E-002 replica against")
    ours = sample_e002()
    for fn in files:
        seen, roots = set(), []
        for line in open(fn):
            if line.startswith("#"):
                continue
            r = int(line.split()[0])
            if r not in seen:
                seen.add(r)
                roots.append(r)
        if roots != ours:
            bad = next((i for i, (a, b) in enumerate(zip(roots, ours)) if a != b), None)
            raise SystemExit(
                f"ERROR: the E-002 replica does not match {os.path.basename(fn)}"
                + (f" at index {bad}: file={roots[bad]} replica={ours[bad]}"
                   if bad is not None else f" (lengths {len(roots)} vs {len(ours)})"))
    if not quiet:
        print(f"E-002 replica matches the root column of all {len(files)} "
              f"committed data/q5_*.txt files")
    return len(files)


def verify_against_c():
    """compile the C sampler's own constants and compare, root by root"""
    with tempfile.TemporaryDirectory() as d:
        src, exe = os.path.join(d, "p.c"), os.path.join(d, "p")
        open(src, "w").write(C_PROBE)
        subprocess.run(["gcc", "-O2", "-o", exe, src], check=True)
        got = [int(x) for x in subprocess.run([exe], capture_output=True,
                                              text=True, check=True).stdout.split()]
    ours = sample_e002()
    if got != ours:
        first = next(i for i, (a, b) in enumerate(zip(got, ours)) if a != b)
        raise SystemExit(f"ERROR: splitmix64 replica diverges from the C at index "
                         f"{first}: C={got[first]} python={ours[first]}")
    print("splitmix64 replica matches the C sampler on all 300 roots\n")


def counts_at(roots, cutoff):
    """how many of `roots` reach a cycle within `cutoff` steps"""
    n = 0
    for u in roots:
        v = u
        for _ in range(cutoff):
            if v in CYCLE5:
                n += 1
                break
            v = 5 * v + 1
            while v % 2 == 0:
                v //= 2
    return n


def orbit(u, steps=STEPS):
    """(steps_to_cycle or None, value_after)"""
    for k in range(steps):
        if u in CYCLE5:
            return k, u
        u = 5 * u + 1
        while u % 2 == 0:
            u //= 2
    # entering at exactly the cutoff is step `steps`, not step 0; reporting 0
    # corrupted the printed "deepest entry into a cycle" (R20-09)
    return (steps if u in CYCLE5 else None), u


def report(label, roots, backs):
    cyc, dig, worst = [], [], 0
    for u in roots:
        k, val = orbit(u)
        if k is None:
            dig.append(len(str(val)))
        else:
            cyc.append(u)
            worst = max(worst, k)
    print(f"{label}  ({backs})")
    print(f"  reach a known cycle : {len(cyc)}")
    print(f"  still growing       : {len(dig)}")
    print(f"  digit range         : {min(dig)} to {max(dig)}")
    # This check has now been wrong twice, in two different ways, and the second
    # way is the instructive one.
    #
    # R18-03: `worst < STEPS` was a tautology, since orbit() reports None past the
    # cutoff; it held at every cutoff down to 1, while the count it guarded moved
    # 19 -> 0.
    #
    # R19-01: the replacement recounted at `2*worst`, which is forced whenever
    # `2*worst <= STEPS`, because `worst` is derived from the very run being
    # audited. Worse, it probed BELOW the cutoff. The only way a cutoff can
    # corrupt this count is by hiding a root that enters a cycle ABOVE it, and no
    # probe below can see that. Constructed input made the old check print
    # "so 5000 is not binding" in a case where 5000 was exactly what bound.
    #
    # So probe ABOVE the cutoff, at a bound fixed independently of the data. This
    # fails exactly when a root enters a cycle between STEPS and the probe, which
    # is the failure the sentence in the paper excludes.
    deep = 2 * STEPS
    at_deep = counts_at(roots, deep)
    assert at_deep == len(cyc), (
        f"the {STEPS}-step cutoff IS binding: {len(cyc)} roots reach a cycle "
        f"within it, but {at_deep} do within {deep}")
    print(f"  deepest entry into a cycle: {worst} steps; extending the cutoff to "
          f"{deep} finds no further root, so {STEPS} is not binding")
    return set(roots), len(cyc)


def main():
    # BOTH checks run on the declared command, not behind a flag (R18-01)
    verify_against_data()
    if "--verify-c" in sys.argv:
        verify_against_c()
    from importlib import import_module
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    "..", "E-001-qx1-empirical-gate"))
    real_e001 = import_module("experiment_gate_richardson").sample_roots(
        random.Random(SEED), N_ROOTS)
    if real_e001 != sample_e001():
        raise SystemExit("ERROR: the E-001 replica no longer matches E-001's own sample_roots")
    print("E-001 replica matches E-001's own sample_roots")

    print(f"accelerated 5x+1, {STEPS} steps per root, "
          f"known cycles {sorted(CYCLE5)}\n")
    a, _ = report("E-001 sample, Python random.Random(2026)",
                  sample_e001(), "backs sections 4 and 5")
    print()
    b, _ = report("E-002 sample, tree_counts.c splitmix64 seed 2026",
                  sample_e002(), "backs section 6 and the data/q5_* files")
    print(f"\nroots common to both samples: {len(a & b)} of {N_ROOTS}")
    print("The two samples are different sets. Section 6's sentence takes the "
          "E-002 numbers.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
