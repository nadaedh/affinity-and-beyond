"""Turn AbMPNN's 100 designed sequences into per-substitution scores.

ProteinMPNN emits whole sequences, not mutation scores. We convert the
sampling frequency at each designed position into a log-odds score against
the parent residue:

    score = log2( (n_mut + 0.5) / (n_wt + 0.5) )

positive  -> AbMPNN prefers the substitution to the parent residue
Fixed positions (everything outside zones A/B) come back as the parent in
every sample and are dropped.

Writes results/abmpnn_scan.csv
"""
import csv
import math
import os
import sys
from collections import Counter

FA = "results/abmpnn/seqs/fv_her2_domainIV.fa"
VH = "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
VL = "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
VH_LEN, LINKER_LEN = 120, 15

ZONE_A = {33, 50, 52, 57, 58, 59, 99, 101, 102, 103, 104, 105,
          163, 165, 166, 167, 184, 185, 188, 201, 226, 227, 228, 229}
ZONE_B = {29, 30, 31, 32, 34, 35, 47, 48, 49, 51, 53, 54, 55, 56, 60, 61, 70,
          97, 98, 100, 106, 107, 108,
          137, 161, 162, 164, 168, 169, 181, 182, 183, 186, 187, 189, 190,
          203, 204, 206, 224, 225, 230, 231}

records = []
with open(FA) as fh:
    head, seq = None, []
    for line in fh:
        line = line.strip()
        if line.startswith(">"):
            if head is not None:
                records.append((head, "".join(seq)))
            head, seq = line, []
        else:
            seq.append(line)
    records.append((head, "".join(seq)))

native_head, native = records[0]
designs = [s for h, s in records[1:]]
print(f"{len(designs)} designed sequences, 1 native reference")

# chain order in the output is alphabetical: A (VL) then B (VH)
nat_vl, nat_vh = native.split("/")
if nat_vl != VL or nat_vh != VH:
    sys.exit("native sequence in the FASTA does not match the parent - check chain order")
print("native reference matches parent VL and VH (chain order A=VL, B=VH confirmed)")

counts = {}          # scfv_pos -> Counter of sampled residues
for d in designs:
    vl, vh = d.split("/")
    for local, aa in enumerate(vh, 1):
        counts.setdefault(local, Counter())[aa] += 1
    for local, aa in enumerate(vl, 1):
        counts.setdefault(local + VH_LEN + LINKER_LEN, Counter())[aa] += 1

parent = {}
for local, aa in enumerate(VH, 1):
    parent[local] = aa
for local, aa in enumerate(VL, 1):
    parent[local + VH_LEN + LINKER_LEN] = aa

rows = []
designed = ZONE_A | ZONE_B
for p in sorted(counts):
    if p not in designed:
        continue
    c, wt = counts[p], parent[p]
    n_wt = c.get(wt, 0)
    for mut, n in c.items():
        if mut == wt:
            continue
        rows.append(dict(scfv_pos=p, wt=wt, mut=mut,
                         zone="A" if p in ZONE_A else "B",
                         n=n, recovery=round(n_wt / len(designs), 3),
                         abmpnn=round(math.log2((n + 0.5) / (n_wt + 0.5)), 4)))

rows.sort(key=lambda r: -r["abmpnn"])
os.makedirs("results", exist_ok=True)
with open("results/abmpnn_scan.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["scfv_pos", "wt", "mut", "zone", "n", "recovery", "abmpnn"])
    w.writeheader()
    w.writerows(rows)

print(f"wrote results/abmpnn_scan.csv ({len(rows)} substitutions seen at "
      f"{len({r['scfv_pos'] for r in rows})} designed positions)")

pos = [r for r in rows if r["abmpnn"] > 0]
print(f"\nAbMPNN prefers the substitution over parent at {len(pos)} (pos,mut) pairs")
print(f"  zone A: {sum(r['zone']=='A' for r in pos)}   zone B: {sum(r['zone']=='B' for r in pos)}")

print(f"\ntop 20 (n = times sampled out of {len(designs)}):")
print(f"{'pos':>5} {'mut':9} {'zone':5} {'n':>4} {'wt kept':>8} {'score':>7}")
for r in rows[:20]:
    print(f"{r['scfv_pos']:5} {r['wt']}->{r['mut']:<6} {r['zone']:5} {r['n']:4} "
          f"{r['recovery']:8.2f} {r['abmpnn']:7.2f}")

low = sorted({(r["scfv_pos"], r["wt"], r["recovery"]) for r in rows}, key=lambda x: x[2])[:12]
print("\nleast-conserved designed positions (AbMPNN most willing to change):")
for p, wt, rec in low:
    print(f"  {wt}{p}: parent kept in {rec:.0%} of samples")
