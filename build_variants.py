#!/usr/bin/env python3
"""Build the 15-variant submission pool for the Adaptyv HER2 scFv challenge."""
import json, csv

d = json.load(open("parent.json"))
PARENT = d["full_sequence"]
DNM = set(d["do_not_mutate"]["positions_scfv"])
LINKER = range(121, 136)

# mutation blocks (scFv 1-based coords)
BASE = ["D62E", "N84S"]                      # deamidation + isomerisation, distal, FR3
STAR = "D102S"                               # DG isomerisation fix, 5.08 A, both models +

VARIANTS = {
 "V01_safe_base":        BASE + ["A72R", "K30S", "S97A"],
 "V02_base_D102S":       BASE + ["A72R", "K30S", "S97A", STAR],
 "V03_safe_alt":         BASE + ["T54D", "A79V", "S97A"],
 "V04_ss_V168L":         BASE + [STAR, "V168L", "A72R", "K30S"],
 "V05_ss_A106Y":         BASE + [STAR, "A106Y", "A72R", "K30S"],
 "V06_ss_combined":      BASE + [STAR, "V168L", "A106Y", "A72R", "K30S"],
 "V07_deamid_ss":        BASE + [STAR, "N55S", "V164I", "A72R"],
 "V08_ox_M107F":         BASE + [STAR, "M107F", "A72R", "K30S"],
 "V09_ox_M107L":         BASE + [STAR, "M107L", "A72R", "K30S"],
 "V10_interface_stack":  BASE + [STAR, "M107F", "A169S", "V168L", "A72R"],
 "V11_charge_D102R":     BASE + ["D102R", "A72R", "K30S", "S97A"],
 "V12_charge_double":    BASE + ["D102R", "G101R", "A72R", "K30S", "S97A"],
 "V13_germline_test":    BASE + ["F188S", "A72R", "K30S", "S97A"],
 "V14_paratope_bet":     BASE + [STAR, "A167Y", "A72R", "K30S"],
 "V15_combined_best":    BASE + [STAR, "M107F", "V168L", "A106Y", "A72R", "K30S"],
}

def apply_muts(seq, muts):
    s = list(seq)
    for m in muts:
        wt, pos, mu = m[0], int(m[1:-1]), m[-1]
        assert s[pos-1] == wt, f"{m}: parent has {s[pos-1]} at {pos}, not {wt}"
        assert pos not in DNM, f"{m}: position {pos} is on the do-not-mutate list"
        assert pos not in LINKER, f"{m}: position {pos} is in the linker"
        s[pos-1] = mu
    return "".join(s)

rows, fasta = [], []
for name, muts in VARIANTS.items():
    assert len(muts) == len(set(int(m[1:-1]) for m in muts)), f"{name}: duplicate position"
    seq = apply_muts(PARENT, muts)
    n = sum(a != b for a, b in zip(PARENT, seq))
    assert n == len(muts), f"{name}: {n} diffs vs {len(muts)} muts"
    assert n >= 5, f"{name}: only {n} mutations, minimum is 5"
    assert len(seq) == 242 and len(seq) <= 260
    assert seq[120:135] == "GGGGSGGGGSGGGGS", f"{name}: linker altered"
    rows.append({"variant": name, "n_mutations": n,
                 "mutations": " ".join(sorted(muts, key=lambda m: int(m[1:-1]))),
                 "length": len(seq), "sequence": seq})
    fasta.append(f">{name} | {n} mut | {' '.join(sorted(muts, key=lambda m: int(m[1:-1])))}\n{seq}")

with open("her2_scfv_15_variants.fasta", "w") as f:
    f.write("\n".join(fasta) + "\n")
with open("her2_scfv_15_variants.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader(); w.writerows(rows)

print(f"{len(rows)} variants built, all checks passed")
print(f"{'variant':24} {'n':>2}  mutations")
for r in rows:
    print(f"{r['variant']:24} {r['n_mutations']:>2}  {r['mutations']}")

# pairwise distinctness
seqs = [r["sequence"] for r in rows]
assert len(set(seqs)) == 15, "duplicate sequences"
print("\nAll 15 sequences distinct. Parent unchanged in pool:", PARENT not in seqs)
