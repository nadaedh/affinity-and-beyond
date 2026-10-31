#!/usr/bin/env python3
"""
Build the truncated HER2 target from 1N8Z for cofolding.

Usage:  python make_target.py [start] [end]       (default 500 640)

Does three things the agent kept stalling on:
  1. finds every chain C residue contacting chains A/B
  2. finds all disulfides in chain C and EXTENDS the window until every
     cysteine inside has its partner inside
  3. writes her2_target.fasta + her2_target.pdb and reports gaps

No external deps beyond biopython.   pip install biopython
"""
import sys, urllib.request, os
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Polypeptide import protein_letters_3to1 as t31

START = int(sys.argv[1]) if len(sys.argv) > 1 else 500
END   = int(sys.argv[2]) if len(sys.argv) > 2 else 640

if not os.path.exists("1n8z.pdb"):
    print("downloading 1N8Z ...")
    urllib.request.urlretrieve("https://files.rcsb.org/download/1N8Z.pdb", "1n8z.pdb")

model = PDBParser(QUIET=True).get_structure("x", "1n8z.pdb")[0]
C = [r for r in model["C"] if r.id[0] == " "]
AB = [r for c in ("A", "B") for r in model[c] if r.id[0] == " "]
cnum = {r.id[1]: r for r in C}

# ---------- 1. contacts ----------
ab_xyz = np.array([a.coord for r in AB for a in r if a.element != "H"])
contacts = []
for r in C:
    xyz = np.array([a.coord for a in r if a.element != "H"])
    d = np.linalg.norm(xyz[:, None] - ab_xyz[None, :], axis=-1).min()
    if d <= 5.0:
        contacts.append((r.id[1], t31.get(r.get_resname(), "X"), round(float(d), 2)))

print(f"\n--- contacts (<=5.0 A to chains A/B): {len(contacts)} residues ---")
print("  " + ", ".join(f"{aa}{n}" for n, aa, _ in contacts))
lo_c, hi_c = contacts[0][0], contacts[-1][0]
print(f"  epitope spans {lo_c}-{hi_c}")

# ---------- 2. disulfides, then extend window to keep pairs whole ----------
cys = [r for r in C if r.get_resname() == "CYS" and "SG" in r]
pairs = []
for i, a in enumerate(cys):
    for b in cys[i + 1:]:
        if np.linalg.norm(a["SG"].coord - b["SG"].coord) < 2.5:
            pairs.append((a.id[1], b.id[1]))
print(f"\n--- disulfides in chain C: {len(pairs)} ---")
print("  " + ", ".join(f"{i}-{j}" for i, j in pairs))

start, end = START, END
for _ in range(50):
    broken = [(i, j) for i, j in pairs
              if (start <= i <= end) != (start <= j <= end)]
    if not broken:
        break
    for i, j in broken:
        start, end = min(start, i, j), max(end, i, j)
    print(f"  extending to {start}-{end} to keep {broken} whole")
else:
    print("  WARNING: could not resolve disulfide pairing")

inside = [(i, j) for i, j in pairs if start <= i <= end and start <= j <= end]
print(f"  final window {start}-{end}; {len(inside)} intact disulfides, 0 orphaned")
assert lo_c >= start and hi_c <= end, "epitope falls outside the window"
print(f"  epitope {lo_c}-{hi_c} sits inside with "
      f"{lo_c-start} / {end-hi_c} residues of margin")

# ---------- 3. write ----------
present = sorted(n for n in cnum if start <= n <= end)
gaps = [(a, b) for a, b in zip(present, present[1:]) if b - a > 1]
if gaps:
    print(f"\n--- unresolved loops in the crystal (modelled by Boltz anyway) ---")
    for a, b in gaps:
        print(f"  gap {a+1}-{b-1}  ({b-a-1} residues)")

seq = "".join(t31.get(cnum[n].get_resname(), "X") for n in present)
with open("her2_target.fasta", "w") as f:
    f.write(f">HER2_domainIV_{start}_{end} | {len(seq)} aa | from 1N8Z chain C\n{seq}\n")

class Sel(Select):
    def accept_chain(self, c):   return c.id == "C"
    def accept_residue(self, r): return r.id[0] == " " and start <= r.id[1] <= end
io = PDBIO(); io.set_structure(model); io.save("her2_target.pdb", Sel())

print(f"\nwrote her2_target.fasta ({len(seq)} aa) and her2_target.pdb")
print(f"complex size with the 242 aa scFv: {len(seq)+242} residues "
      f"(vs {len(C)+242} untruncated)")
