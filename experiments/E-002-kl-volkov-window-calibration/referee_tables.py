#!/usr/bin/env python3
"""
E-002, part 9: the two tables the paper prints, generated rather than typed.

An external referee for Experimental Mathematics asked for the argument to be
carried by two tables: one placing every construction against its pressure
target at the calibration depth, and one showing the comparison across every
checkpoint decade the grid supports. This produces both, in LaTeX and in plain
text, from the committed count files.

The second calibrator at 0.678 (`iidq`) exists because of that referee: until
this round the 0.678 target had a single construction, so the calibration band
measured on the 0.650919 family was being carried to a target where its width
had never been checked. It now has been.

Estimator: `summary.py`'s own `decade_estimate`, imported. The four readings
the paper already published must come out unchanged; the script checks that
and aborts if not.

Serial (Rule 9c): arithmetic over committed count files, about 30 s.

Run: python3 referee_tables.py
"""
import sys

sys.path.insert(0, ".")
from summary import decade_estimate                      # noqa: E402

DEC = 10
# name, pressure target, file stem, structural match to the arithmetic tree
ROWS = [
    ("arithmetic tree",        None,       "q5_arith_b15",   "exact object"),
    ("integer recursion",      0.650919,   "q5_cyc_b15",     "identity and congruence"),
    ("real-valued relaxation", 0.650919,   "q5_cycq500_b15", "identity up to constant, congruence"),
    ("unconstrained recursion", 0.650919,  "q5_iid_b15",     "identity, no congruence"),
    ("unconstrained relaxation", 0.650919, "q5_iidq500_b15", "identity up to constant, no congruence"),
    ("relaxation retuned",     0.678,      "q5_cycq505_b15", "identity up to constant, congruence"),
    ("unconstrained retuned",  0.678,      "q5_iidq505_b15", "identity up to constant, no congruence"),
]
PUBLISHED = {"q5_arith_b15": 0.64926, "q5_cyc_b15": 0.65122,
             "q5_cycq500_b15": 0.64981, "q5_cycq505_b15": 0.67748}


def main():
    v = {}
    for _, _, stem, _ in ROWS:
        p, lo, hi = decade_estimate(f"data/{stem}.txt", dec_top=DEC)
        v[stem] = (p, lo, hi)
        if stem in PUBLISHED and round(p, 5) != PUBLISHED[stem]:
            raise SystemExit(f"ERROR: {stem} reads {p:.5f}, the paper publishes "
                             f"{PUBLISHED[stem]:.5f}")
    a = v["q5_arith_b15"][0]

    print(f"TABLE 1: every construction at the calibration decade "
          f"1e{DEC-1} -> 1e{DEC}, grid b15\n")
    print(f"{'construction':<28}{'target':>10}{'reading':>10}{'offset':>10}"
          f"{'95% interval':>20}  structural match")
    for name, tgt, stem, match in ROWS:
        p, lo, hi = v[stem]
        off = "" if tgt is None else f"{p - a:+.5f}"
        tg = "unknown" if tgt is None else f"{tgt:.6f}"
        print(f"{name:<28}{tg:>10}{p:>10.5f}{off:>10}"
              f"   [{lo:.5f},{hi:.5f}]  {match}")

    for tgt in (0.650919, 0.678):
        band = [v[s][0] for _, t, s, _ in ROWS if t == tgt]
        n = len(band)
        print(f"\n  target {tgt}: {n} constructions, band [{min(band):.5f}, "
              f"{max(band):.5f}], width {max(band)-min(band):.5f}")
        if tgt == 0.650919:
            inside = min(band) <= a <= max(band)
            print(f"    the arithmetic tree reads {a:.5f}, "
                  f"{'INSIDE' if inside else 'OUTSIDE'} this band")
        else:
            print(f"    the arithmetic tree reads {a:.5f}, "
                  f"{min(band)-a:.5f} below the floor of this band")

    print(f"\n\nTABLE 2: the comparison at every checkpoint decade the grid supports\n")
    print(f"{'decade':<12}{'arithmetic':>12}{'0.650919 band':>26}"
          f"{'0.678 band':>26}{'separation':>12}")
    for dec in (6, 7, 8, 9, 10):
        r = {}
        for _, _, stem, _ in ROWS:
            r[stem], _, _ = decade_estimate(f"data/{stem}.txt", dec_top=dec)
        aa = r["q5_arith_b15"]
        b0 = [r[s] for _, t, s, _ in ROWS if t == 0.650919]
        b7 = [r[s] for _, t, s, _ in ROWS if t == 0.678]
        sep = min(b7) - aa
        print(f"1e{dec-1}->1e{dec:<7}{aa:>12.5f}"
              f"   [{min(b0):.5f},{max(b0):.5f}]"
              f"   [{min(b7):.5f},{max(b7):.5f}]{sep:>12.5f}")
    print("\n  separation is measured to the FLOOR of the 0.678 band, which is the")
    print("  conservative direction: it is the smallest gap the two families admit.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
