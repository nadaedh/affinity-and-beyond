"""Interface residues of the trastuzumab Fab against HER2 in 1N8Z.

Restricted to the variable domains, since only those are in the scFv.
Two definitions reported together:
  contact  -- heavy atom within CUTOFF of a HER2 heavy atom
  dSASA    -- SASA lost by that residue when HER2 is added (A^2)
"""
import copy
import sys
from Bio.PDB import MMCIFParser, NeighborSearch
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1

CUTOFF = 4.5
VH_LEN, LINKER_LEN, VL_LEN = 120, 15, 107
PARENT = ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
          "GGGGSGGGGSGGGGS"
          "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK")

# chain -> (role, last residue number that is part of the scFv)
CHAINS = {"B": ("VH", VH_LEN), "A": ("VL", VL_LEN)}

cif = sys.argv[1]
model = next(MMCIFParser(QUIET=True).get_structure("x", cif).get_models())
her2 = [r for r in model["C"] if is_aa(r, standard=True)]
ns = NeighborSearch([a for r in her2 for a in r if a.element != "H"])
sr = ShrakeRupley()


def scfv_index(chain, num):
    if chain == "B":
        return num
    return num + VH_LEN + LINKER_LEN


rows = []
for ch, (role, last) in CHAINS.items():
    residues = [r for r in model[ch] if is_aa(r, standard=True) and r.id[1] <= last]

    # SASA of this variable domain alone, and with HER2 present
    alone = copy.deepcopy(model)
    for cid in [c.id for c in alone]:
        if cid != ch:
            alone.detach_child(cid)
    for r in list(alone[ch]):
        if r.id[1] > last or not is_aa(r, standard=True):
            alone[ch].detach_child(r.id)
    sr.compute(alone, level="R")
    free = {r.id[1]: r.sasa for r in alone[ch]}

    withC = copy.deepcopy(alone)
    withC.add(copy.deepcopy(model["C"]))
    sr.compute(withC, level="R")
    bound = {r.id[1]: r.sasa for r in withC[ch]}

    for r in residues:
        num, aa = r.id[1], seq1(r.get_resname())
        partners = set()
        for a in r:
            if a.element != "H":
                for hit in ns.search(a.coord, CUTOFF):
                    partners.add(hit.get_parent().id[1])
        d = free.get(num, 0) - bound.get(num, 0)
        if partners or d >= 1.0:
            rows.append((role, num, aa, sorted(partners), d))

print(f"Trastuzumab paratope vs HER2 (1N8Z), contact cutoff {CUTOFF} A\n")
print(f"{'dom':4} {'res':6} {'scFv':>5} {'dSASA':>7}  HER2 contacts")
print("-" * 64)
mismatch = 0
for role, num, aa, partners, d in sorted(rows, key=lambda x: (x[0] != "VH", x[1])):
    i = scfv_index("B" if role == "VH" else "A", num)
    ok = PARENT[i - 1] == aa
    mismatch += not ok
    flag = "" if ok else f" !!{PARENT[i-1]}"
    p = ",".join(str(x) for x in partners[:7]) or "-"
    print(f"{role:4} {aa}{num:<5} {i:5d} {d:7.1f}  {p}{flag}")

print(f"\n{len(rows)} interface residues; sequence mismatches vs parent: {mismatch}")
contact_only = [r for r in rows if r[3]]
print(f"{len(contact_only)} make direct contact (<{CUTOFF} A)")
