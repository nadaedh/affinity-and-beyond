#!/usr/bin/env python3
"""
Score every Boltz-2 prediction in pred/ : ipSAE + binding-mode RMSD vs 1N8Z.

  python score_all.py

Needs: ipsae (pip install ipsae), biopython, 1n8z.pdb, pred/
Writes: results.csv  and prints a ranked table.

RMSD = superpose the predicted target chain onto 1N8Z chain C, then measure
CA deviation of the scFv. Low RMSD means the designed binding mode matches
the crystal.
"""
import glob, os, re, subprocess, csv, sys
import numpy as np
from Bio.PDB import PDBParser, MMCIFParser, Superimposer

TARGET_START = 489          # predicted chain B residue 1 == HER2 489
VH_LEN, LINK_LEN = 120, 15

ref = PDBParser(QUIET=True).get_structure("r", "1n8z.pdb")[0]
ref_C  = {r.id[1]: r for r in ref["C"] if r.id[0] == " " and "CA" in r}
ref_VH = {r.id[1]: r for r in ref["B"] if r.id[0] == " " and "CA" in r and r.id[1] <= VH_LEN}
ref_VL = {r.id[1]: r for r in ref["A"] if r.id[0] == " " and "CA" in r and r.id[1] <= 107}

def chains_of(cif):
    m = MMCIFParser(QUIET=True).get_structure("m", cif)[0]
    ch = [c for c in m]
    ch.sort(key=lambda c: -len([r for r in c if r.id[0] == " "]))
    return ch[0], ch[1]           # longer = scFv (242), shorter = target (119)

def rmsd_vs_crystal(cif):
    try:
        scfv, tgt = chains_of(cif)
        t = [r for r in tgt if r.id[0] == " " and "CA" in r]
        fixed, moving = [], []
        for i, r in enumerate(t):                       # superpose on target
            n = TARGET_START + i
            if n in ref_C:
                fixed.append(ref_C[n]["CA"]); moving.append(r["CA"])
        if len(fixed) < 30: return None, None
        sup = Superimposer(); sup.set_atoms(fixed, moving)
        tgt_rmsd = round(sup.rms, 2)

        s = [r for r in scfv if r.id[0] == " " and "CA" in r]
        sup.apply([a for r in s for a in r])             # move scFv with it
        d = []
        for i, r in enumerate(s, start=1):
            if i <= VH_LEN and i in ref_VH:
                d.append(r["CA"].coord - ref_VH[i]["CA"].coord)
            elif i > VH_LEN + LINK_LEN:
                j = i - VH_LEN - LINK_LEN
                if j in ref_VL:
                    d.append(r["CA"].coord - ref_VL[j]["CA"].coord)
        if not d: return tgt_rmsd, None
        return tgt_rmsd, round(float(np.sqrt((np.array(d) ** 2).sum(1).mean())), 2)
    except Exception as e:
        print("  rmsd failed:", e); return None, None

def run_ipsae(npz, cif):
    try:
        subprocess.run(["ipsae", npz, cif, "10", "10"],
                       capture_output=True, text=True, timeout=300)
    except Exception as e:
        print("  ipsae failed:", e); return None
    base = cif[:-4]
    for pat in (base + "*_chains.csv", os.path.dirname(cif) + "/*chains*.csv",
                os.path.dirname(cif) + "/*chains*.txt"):
        for f in glob.glob(pat):
            txt = open(f).read()
            vals = [float(x) for x in re.findall(r"0\.\d+", txt)]
            if vals: return round(max(vals), 4)
    return None

rows = []
for d in sorted(glob.glob("pred/*/")):
    name = os.path.basename(d.rstrip("/"))
    cifs = glob.glob(d + "**/predictions/**/*_model_0.cif", recursive=True)
    npzs = glob.glob(d + "**/predictions/**/pae_*_model_0.npz", recursive=True)
    if not cifs or not npzs:
        print(f"{name:22} no output yet"); continue
    print(f"{name:22} scoring ...")
    ipsae = run_ipsae(npzs[0], cifs[0])
    tgt_rmsd, scfv_rmsd = rmsd_vs_crystal(cifs[0])
    rows.append({"variant": name, "ipSAE": ipsae,
                 "target_fit_rmsd": tgt_rmsd, "scFv_rmsd_vs_1N8Z": scfv_rmsd})

if not rows:
    sys.exit("no predictions scored")

base = next((r["ipSAE"] for r in rows if r["variant"] == "PARENT"), None)
for r in rows:
    r["delta_vs_parent"] = (round(r["ipSAE"] - base, 4)
                            if base and r["ipSAE"] is not None else None)

with open("results.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

rows.sort(key=lambda r: (r["ipSAE"] is None, -(r["ipSAE"] or 0)))
print(f"\n{'variant':22} {'ipSAE':>7} {'d_parent':>9} {'scFv RMSD':>10}")
print("-" * 52)
for r in rows:
    print(f"{r['variant']:22} {str(r['ipSAE']):>7} "
          f"{str(r['delta_vs_parent']):>9} {str(r['scFv_rmsd_vs_1N8Z']):>10}")
print(f"\nparent baseline ipSAE: {base}")
print("wrote results.csv")
