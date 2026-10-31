"""Zone assignment (A/B/C) and developability liability motifs for the parent scFv.

Zone A: direct contact with HER2 (heavy atom < 4.5 A)
Zone B: second shell -- within 5 A of a zone-A side chain, but not zone A itself
Zone C: everything else
"""
import copy
import re
import sys
from Bio.PDB import MMCIFParser, NeighborSearch
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1

CONTACT, SHELL = 4.5, 5.0
VH_LEN, LINKER_LEN, VL_LEN = 120, 15, 107
PARENT = ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
          "GGGGSGGGGSGGGGS"
          "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK")
CHAINS = {"B": ("VH", VH_LEN, 0), "A": ("VL", VL_LEN, VH_LEN + LINKER_LEN)}

model = next(MMCIFParser(QUIET=True).get_structure("x", sys.argv[1]).get_models())
her2_atoms = [a for r in model["C"] if is_aa(r, standard=True) for a in r if a.element != "H"]
ns_her2 = NeighborSearch(her2_atoms)

# collect variable-domain residues, keyed by scFv index
vres = {}
for ch, (role, last, off) in CHAINS.items():
    for r in model[ch]:
        if is_aa(r, standard=True) and r.id[1] <= last:
            vres[off + r.id[1]] = (role, r)

zoneA = set()
for i, (role, r) in vres.items():
    for a in r:
        if a.element != "H" and ns_her2.search(a.coord, CONTACT):
            zoneA.add(i)
            break

# zone B: near a zone-A side chain but not in contact itself
a_atoms = [a for i in zoneA for a in vres[i][1] if a.element != "H"]
ns_a = NeighborSearch(a_atoms)
zoneB = set()
for i, (role, r) in vres.items():
    if i in zoneA:
        continue
    for a in r:
        if a.element != "H" and ns_a.search(a.coord, SHELL):
            zoneB.add(i)
            break


def zone(i):
    return "A" if i in zoneA else ("B" if i in zoneB else "C")


print(f"Zone A (contact): {len(zoneA)}   Zone B (second shell): {len(zoneB)}   "
      f"Zone C: {len(vres) - len(zoneA) - len(zoneB)}  of {len(vres)} variable-domain residues\n")

print("Zone B positions (scFv numbering):")
for role in ("VH", "VL"):
    pos = sorted(i for i in zoneB if vres[i][0] == role)
    print(f"  {role}: " + " ".join(f"{PARENT[i-1]}{i}" for i in pos))

# --- developability liability motifs ---
MOTIFS = [("NG", "deamidation, fast"), ("NS", "deamidation"), ("NT", "deamidation, slow"),
          ("DG", "isomerisation, fast"), ("DS", "isomerisation"), ("DD", "isomerisation")]
print("\n--- liability motifs in the 242 aa parent ---")
print(f"{'motif':6} {'pos':>8}  {'zones':8} {'risk':22} burial")
print("-" * 62)
for m, risk in MOTIFS:
    for hit in re.finditer(f"(?={m})", PARENT):
        i = hit.start() + 1          # 1-based position of first residue
        if i not in vres and i + 1 not in vres:
            continue                  # inside the linker
        zs = f"{zone(i)}/{zone(i+1)}"
        flag = "  <-- IN CONTACT" if "A" in zs else ""
        print(f"{m:6} {PARENT[i-1]}{i}-{PARENT[i]}{i+1:<4} {zs:8} {risk:22}{flag}")

print("\nfree cysteines:", [f"C{i}" for i, c in enumerate(PARENT, 1) if c == "C"])
