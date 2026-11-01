"""Build the trimmed Fv-HER2 complex used for AbMPNN design and Boltz cofolding.

Keeps VH (chain B 1-120), VL (chain A 1-107) and HER2 domain IV (chain C),
trimmed around the observed epitope with a margin.

Outputs:
  data/structures/fv_her2_domainIV.pdb   - for AbMPNN
  data/her2_domainIV.fasta               - target sequence for Boltz
  results/abmpnn_design_positions.txt    - zone A+B positions, ProteinMPNN format
"""
import os

from Bio.PDB import MMCIFParser, PDBIO, Select
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1

VH_LAST, VL_LAST = 120, 107
HER2_LO, HER2_HI = 490, 620          # domain IV; observed contacts are 557-605
VH_LEN, LINKER_LEN = 120, 15

ZONE_A = {33, 50, 52, 57, 58, 59, 99, 101, 102, 103, 104, 105,
          163, 165, 166, 167, 184, 185, 188, 201, 226, 227, 228, 229}
ZONE_B = {29, 30, 31, 32, 34, 35, 47, 48, 49, 51, 53, 54, 55, 56, 60, 61, 70,
          97, 98, 100, 106, 107, 108,
          137, 161, 162, 164, 168, 169, 181, 182, 183, 186, 187, 189, 190,
          203, 204, 206, 224, 225, 230, 231}
FORBIDDEN = {22, 96, 158, 223}       # conserved cysteines - never design


class Trim(Select):
    def accept_residue(self, r):
        if not is_aa(r, standard=True):
            return 0
        ch, n = r.get_parent().id, r.id[1]
        if ch == "B":
            return n <= VH_LAST
        if ch == "A":
            return n <= VL_LAST
        if ch == "C":
            return HER2_LO <= n <= HER2_HI
        return 0


model = next(MMCIFParser(QUIET=True).get_structure("x", "data/structures/1n8z.cif").get_models())
for cid in [c.id for c in model]:
    if cid not in ("A", "B", "C"):
        model.detach_child(cid)

os.makedirs("results", exist_ok=True)
io = PDBIO()
io.set_structure(model)
out_pdb = "data/structures/fv_her2_domainIV.pdb"
io.save(out_pdb, select=Trim())

kept = {c.id: [r for r in model[c.id] if Trim().accept_residue(r)] for c in model}
for cid in ("B", "A", "C"):
    res = kept[cid]
    print(f"chain {cid}: {len(res)} residues, {res[0].id[1]}-{res[-1].id[1]}")

her2 = "".join(seq1(r.get_resname()) for r in kept["C"])
with open("data/her2_domainIV.fasta", "w") as fh:
    fh.write(f">HER2_domainIV_{HER2_LO}-{HER2_HI}\n{her2}\n")
print(f"\nHER2 domain IV: {len(her2)} aa -> data/her2_domainIV.fasta")

# design positions for AbMPNN, in ProteinMPNN's per-chain local numbering
design = sorted((ZONE_A | ZONE_B) - FORBIDDEN)
vh = [i for i in design if i <= VH_LEN]
vl = [i - VH_LEN - LINKER_LEN for i in design if i > VH_LEN + LINKER_LEN]
with open("results/abmpnn_design_positions.txt", "w") as fh:
    fh.write(f"B: {' '.join(map(str, vh))}\n")
    fh.write(f"A: {' '.join(map(str, vl))}\n")

print(f"design positions: {len(vh)} in VH (chain B), {len(vl)} in VL (chain A)")
print(f"complex written to {out_pdb}")
