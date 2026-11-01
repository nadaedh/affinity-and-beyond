"""Per-position annotation table for the 242 aa scFv.

Locates the six CDRs by exact substring match against the known trastuzumab
loops (no numbering scheme assumed, so nothing to get wrong), assigns zones,
and computes the CB-CB distance matrix used by the Phase 2 combination rules.

Outputs:
  results/positions.csv          one row per position
  results/cb_distances.npy       227x227 CB-CB distances, indexed by scFv pos
"""
import csv
import os

import numpy as np
from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import is_aa

SCFV = ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
        "GGGGSGGGGSGGGGS"
        "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK")
VH_LEN, LINKER_LEN = 120, 15

# canonical trastuzumab CDR sequences (Kabat definition), located by substring
CDRS = {"CDR-H1": "DTYIH", "CDR-H2": "RIYPTNGYTRYADSVKG", "CDR-H3": "WGGDGFYAMDY",
        "CDR-L1": "RASQDVNTAVA", "CDR-L2": "SASFLYS", "CDR-L3": "QQHYTTPPT"}

# Zone A is a DIRECT CONTACT (<4.5 A heavy atom), per the plan's definition.
# 101, 102, 184 and 201 bury surface on binding but make no contact at that
# cutoff, so they are interface RIM -> zone B. An earlier version of these
# scripts wrongly listed them in zone A; positions.csv is the source of truth.
ZONE_A = {33, 50, 52, 57, 58, 59, 99, 103, 104, 105,
          163, 165, 166, 167, 185, 188, 226, 227, 228, 229}
ZONE_B = {29, 30, 31, 32, 34, 35, 47, 48, 49, 51, 53, 54, 55, 56, 60, 61, 70,
          97, 98, 100, 101, 102, 106, 107, 108,
          137, 161, 162, 164, 168, 169, 181, 182, 183, 184, 186, 187, 189, 190,
          201, 203, 204, 206, 224, 225, 230, 231}
FORBIDDEN = {22, 96, 158, 223}

region = {}
for name, seq in CDRS.items():
    n = SCFV.count(seq)
    if n != 1:
        raise SystemExit(f"{name}: expected 1 match for {seq!r}, found {n}")
    start = SCFV.index(seq) + 1
    for i in range(start, start + len(seq)):
        region[i] = name
    print(f"{name:8} {seq:20} scFv {start}-{start + len(seq) - 1}")

# --- CB coordinates from the crystal structure ---
model = next(MMCIFParser(QUIET=True).get_structure("x", "data/structures/1n8z.cif").get_models())
coord = {}
for ch, (last, off) in {"B": (VH_LEN, 0), "A": (107, VH_LEN + LINKER_LEN)}.items():
    for r in model[ch]:
        if not is_aa(r, standard=True) or r.id[1] > last:
            continue
        a = r["CB"] if "CB" in r else r["CA"]      # glycine has no CB
        coord[off + r.id[1]] = a.coord

n = len(SCFV)
D = np.full((n + 1, n + 1), np.nan)
for i, ci in coord.items():
    for j, cj in coord.items():
        D[i, j] = float(np.linalg.norm(ci - cj))

os.makedirs("results", exist_ok=True)
np.save("results/cb_distances.npy", D)

with open("results/positions.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["scfv_pos", "wt", "domain", "region", "zone", "designable"])
    for i, aa in enumerate(SCFV, 1):
        if i <= VH_LEN:
            dom = "VH"
        elif i <= VH_LEN + LINKER_LEN:
            dom = "linker"
        else:
            dom = "VL"
        z = "A" if i in ZONE_A else ("B" if i in ZONE_B else "C")
        reg = region.get(i, "linker" if dom == "linker" else "FR")
        ok = dom != "linker" and i not in FORBIDDEN
        w.writerow([i, aa, dom, reg, z if dom != "linker" else "-", int(ok)])

print(f"\nCB-CB matrix: {len(coord)} positions with coordinates")
print("wrote results/positions.csv and results/cb_distances.npy")

# sanity: zone A residues should mostly sit in CDRs
inA_cdr = sum(1 for i in ZONE_A if region.get(i, "FR").startswith("CDR"))
print(f"zone A in CDRs: {inA_cdr}/{len(ZONE_A)}")
for nm in CDRS:
    za = sorted(i for i in ZONE_A if region.get(i) == nm)
    zb = sorted(i for i in ZONE_B if region.get(i) == nm)
    print(f"  {nm}: zoneA {za}  zoneB {zb}")
