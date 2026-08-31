#!/usr/bin/env python3
"""
E-002 analysis: turn a count matrix produced by `tree_counts` into

  1. the E-001 estimator, reproduced exactly (fixed window, Aitken in the
     truncation buffer, bootstrap over roots);
  2. a per-decade slope panel in which EVERY decade is buffer-extrapolated
     separately before the decade sequence is read, which the E-001 panel
     was not;
  3. the deficit of that per-decade slope against alpha_-(q), decade by
     decade. On a control run alpha_-(q) is the exponent the process
     is TUNED TO, so the deficit column is the distance from that target,
     measured directly. Reading it as estimator bias requires the
     realized-exponent assumption, which the paper states and does not
     prove. On the arithmetic run the same column is a distance and
     nothing more, since the arithmetic exponent is what is in dispute.

No correction-form fit is reported. The exact annealed count
(`annealed_exact.py`) reaches slope 0.6517 already at t=3, so the residual
here is not an annealed correction whose exponent could be fitted; it is the
lag of a single realization's log-slope behind the mean, and neither a
power-law nor a 1/L ansatz was derived for that.

Usage: python3 analyze.py data/<file>.txt [--label NAME] [--truth VALUE]
"""
import math
import random
import sys


def load(path):
    roots, counts = [], []
    hdr = {}
    resid = None
    for line in open(path):
        if line.startswith("#resid"):
            resid = [int(x) for x in line.split()[1:]]
            continue
        if line.startswith("#"):
            for kv in line[1:].split():
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    hdr[k] = v
            continue
        f = line.split()
        roots.append(int(f[0]))
        counts.append([int(x) for x in f[2:]])
    cp_lo, cp_hi = int(hdr["cp_lo"]), int(hdr["cp_hi"])
    buf_lo, buf_hi = int(hdr["buf_lo"]), int(hdr["buf_hi"])
    n_cp, n_buf = cp_hi - cp_lo + 1, buf_hi - buf_lo + 1
    mats = [[row[ci * n_buf:(ci + 1) * n_buf] for ci in range(n_cp)] for row in counts]
    return hdr, roots, mats, cp_lo, buf_lo, n_cp, n_buf, resid


def aitken(s):
    if len(s) < 3:
        return None
    s1, s2, s3 = s[-3], s[-2], s[-1]
    # already converged in the buffer direction: extrapolating a difference
    # that is pure rounding noise would amplify it, so take the last value
    if abs(s3 - s2) < 1e-5:
        return s3
    den = (s3 - s2) - (s2 - s1)
    if abs(den) < 1e-12:
        return None
    return s3 - (s3 - s2) ** 2 / den


def slope(mat, ci_lo, ci_hi, bi):
    lo, hi = mat[ci_lo][bi], mat[ci_hi][bi]
    if lo <= 0 or hi <= 0:
        return None
    return math.log10(hi / lo) / (ci_hi - ci_lo)


def mean(v):
    return sum(v) / len(v)


def alpha_minus(q):
    """smaller root of q^(alpha-1) = 2^alpha - 1"""
    if q == 3:
        return 1.0
    lo, hi = 1e-9, 1.0 - 1e-12
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        if (q ** (mid - 1.0)) - (2.0 ** mid - 1.0) > 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def main():
    path = sys.argv[1]
    label = sys.argv[sys.argv.index("--label") + 1] if "--label" in sys.argv else path
    hdr, roots, mats, cp_lo, buf_lo, n_cp, n_buf, resid = load(path)
    print(f"=== {label} ===")
    print(f"    q={hdr['q']} mode={hdr['mode']} roots={len(roots)} "
          f"cp=1e{cp_lo}..1e{cp_lo+n_cp-1} buf=1e{buf_lo}..1e{buf_lo+n_buf-1}")

    if resid and sum(resid) > 0:   # mode cycq carries no integer residues
        tot = sum(resid)
        print("    residue histogram of counted nodes (u mod q): "
              + " ".join(f"{r}:{c/tot:.5f}" for r, c in enumerate(resid))
              + f"   [uniform = {1/len(resid):.5f}]")

    # ---- 1. the E-001 estimator, on the window 1e5..1e8 if present ----
    if cp_lo <= 5 and cp_lo + n_cp - 1 >= 8:
        a, b = 5 - cp_lo, 8 - cp_lo
        per_root = []
        for m in mats:
            row = [slope(m, a, b, bi) for bi in range(n_buf)]
            if all(v is not None for v in row):
                per_root.append(row)
        curve = [mean([r[bi] for r in per_root]) for bi in range(n_buf)]
        print("\n    [1] E-001 estimator, fixed window 1e5..1e8, mean slope by buffer:")
        print("        " + "  ".join(f"1e{buf_lo+bi}:{curve[bi]:.5f}" for bi in range(n_buf)))
        ait = aitken(curve)
        if ait:
            boots = []
            rng = random.Random(3)
            for _ in range(2000):
                idx = [rng.randrange(len(per_root)) for _ in range(len(per_root))]
                c = [mean([per_root[i][bi] for i in idx]) for bi in range(n_buf)]
                v = aitken(c)
                if v is not None:
                    boots.append(v)
            boots.sort()
            print(f"        Aitken(buffer->inf) = {ait:.5f}   bootstrap 95% CI "
                  f"[{boots[int(0.025*len(boots))]:.5f}, {boots[int(0.975*len(boots))-1]:.5f}]")

    # ---- 2. per-decade slopes, each buffer-extrapolated separately ----
    print("\n    [2] per-decade slope, raw at each buffer, then buffer-extrapolated")
    print("        with a bootstrap over roots (2000 resamples):")
    print("        decade      " + "  ".join(f"1e{buf_lo+bi}" for bi in range(n_buf))
          + "    Aitken   95% CI")
    dec_L, dec_y = [], []
    rng = random.Random(17)
    for ci in range(n_cp - 1):
        # per-root slope at each buffer for this decade
        per_root = []
        for m in mats:
            r = [slope(m, ci, ci + 1, bi) for bi in range(n_buf)]
            if all(v is not None for v in r):
                per_root.append(r)
        if not per_root:
            continue
        row = [mean([r[bi] for r in per_root]) for bi in range(n_buf)]
        # only buffers strictly above the decade's upper checkpoint are usable
        use_idx = [bi for bi in range(n_buf) if buf_lo + bi > cp_lo + ci + 1]
        usable = [row[bi] for bi in use_idx]
        ait = aitken(usable) if len(usable) >= 3 else None
        ci_lo = ci_hi = None
        if ait is not None:
            boots = []
            for _ in range(2000):
                idx = [rng.randrange(len(per_root)) for _ in range(len(per_root))]
                c = [mean([per_root[i][bi] for i in idx]) for bi in use_idx]
                v = aitken(c)
                if v is not None:
                    boots.append(v)
            boots.sort()
            ci_lo = boots[int(0.025 * len(boots))]
            ci_hi = boots[int(0.975 * len(boots)) - 1]
        cells = "  ".join(f"{row[bi]:.4f}" for bi in range(n_buf))
        band = f"[{ci_lo:.4f},{ci_hi:.4f}]" if ci_lo is not None else ""
        print(f"        1e{cp_lo+ci}->1e{cp_lo+ci+1}  {cells}    "
              + (f"{ait:.4f}" if ait is not None else "  --  ")
              + f"  {band}  (n_buf usable {len(usable)}, n_roots {len(per_root)})")
        if ait is not None and len(usable) >= 4:
            dec_L.append(cp_lo + ci + 0.5)
            dec_y.append(ait)

    # ---- 3. deficit against alpha_-(q) ----
    truth = float(sys.argv[sys.argv.index("--truth") + 1]) if "--truth" in sys.argv \
        else alpha_minus(int(hdr["q"]))
    if dec_L:
        kind = "distance from pressure target (bias only under the realized-exponent assumption)" \
            if hdr["mode"] != "arith" else "distance only (the arithmetic exponent is the disputed quantity)"
        print(f"\n    [3] deficit of the buffer-extrapolated decade slope against "
              f"alpha_- = {truth:.6f}")
        print(f"        {kind}")
        print("        " + "  ".join(f"L={L}:{truth-y:+.4f}" for L, y in zip(dec_L, dec_y)))
        print(f"        Kontorovich-Lagarias vs Volkov separation Delta = {0.678-0.650919:.4f}")
    print()


if __name__ == "__main__":
    main()
