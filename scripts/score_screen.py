"""Aggregate the Boltz screen into one ranked table.

For every prediction, reads Boltz's own confidence metrics and computes, from
the predicted structure:

  epitope_jaccard  overlap of the HER2 residues this variant contacts with the
                   ones the PARENT contacts at the same seed. Catches a design
                   that scores well while docking somewhere else.
  pose_ok          epitope_jaccard >= POSE_MIN
  d_iptm           paired delta vs the parent at the SAME seed

If ipSAE output files are present (run scripts/run_ipsae.sh first) their
ipSAE_min is merged in as well.

Writes results/screen_scores.csv
"""
import csv
import glob
import json
import os
import re
from collections import defaultdict

import numpy as np
from Bio.PDB import MMCIFParser, NeighborSearch
from Bio.PDB.Polypeptide import is_aa

CONTACT = 4.5
POSE_MIN = 0.60
ROOT = "results/boltz"
parser = MMCIFParser(QUIET=True)


def epitope(cif):
    """HER2 (chain B) residues contacted by the scFv (chain A)."""
    model = next(parser.get_structure("x", cif).get_models())
    ids = [c.id for c in model]
    if len(ids) < 2:
        return None
    a, b = model[ids[0]], model[ids[1]]
    atoms_a = [at for r in a if is_aa(r, standard=True) for at in r if at.element != "H"]
    if not atoms_a:
        return None
    ns = NeighborSearch(atoms_a)
    hit = set()
    for r in b:
        if not is_aa(r, standard=True):
            continue
        for at in r:
            if at.element != "H" and ns.search(at.coord, CONTACT):
                hit.add(r.id[1])
                break
    return hit


def read_ipsae(d):
    """Parse ipsae.py's <model>_<pae>_<dist>.txt.

    It emits one asym row per chain direction (A->B and B->A). ipSAE_min is
    the smaller of the two, which is the value the thresholds refer to.
    """
    hits = [f for f in glob.glob(os.path.join(d, "*_10_10.txt"))
            if not f.endswith("_byres.txt")]
    for f in hits:
        try:
            lines = [l.split() for l in open(f) if l.strip()]
        except OSError:
            continue
        vals = []
        for parts in lines[1:]:                     # skip header
            if len(parts) > 5 and parts[4] == "asym":
                try:
                    vals.append(float(parts[5]))    # ipSAE column
                except ValueError:
                    pass
        if vals:
            return round(min(vals), 4)              # ipSAE_min
    return None


rows = []
for seeddir in sorted(glob.glob(f"{ROOT}/seed*")):
    seed = int(os.path.basename(seeddir).replace("seed", ""))
    for rundir in sorted(glob.glob(f"{seeddir}/*_s{seed}")):
        name = os.path.basename(rundir)[: -len(f"_s{seed}")]
        pred = glob.glob(f"{rundir}/**/confidence_*_model_0.json", recursive=True)
        if not pred:
            continue
        conf = json.load(open(pred[0]))
        d = os.path.dirname(pred[0])
        cif = glob.glob(f"{d}/*_model_0.cif")
        ep = epitope(cif[0]) if cif else None
        rows.append(dict(name=name, seed=seed,
                         iptm=round(conf["iptm"], 4),
                         ptm=round(conf["ptm"], 4),
                         plddt=round(conf["complex_plddt"], 4),
                         iplddt=round(conf["complex_iplddt"], 4),
                         ipde=round(conf["complex_ipde"], 3),
                         ipsae=read_ipsae(d),
                         epitope=ep))

if not rows:
    raise SystemExit("no finished predictions found under results/boltz/")

# paired comparison against the parent at the same seed
parent = {r["seed"]: r for r in rows if r["name"] == "parent"}
for r in rows:
    p = parent.get(r["seed"])
    r["d_iptm"] = round(r["iptm"] - p["iptm"], 4) if p else None
    if p and r["epitope"] is not None and p["epitope"] is not None:
        u = r["epitope"] | p["epitope"]
        r["epitope_jaccard"] = round(len(r["epitope"] & p["epitope"]) / len(u), 3) if u else 0.0
    else:
        r["epitope_jaccard"] = None
    r["n_epitope"] = len(r["epitope"]) if r["epitope"] is not None else None
    r["pose_ok"] = (r["epitope_jaccard"] is not None and r["epitope_jaccard"] >= POSE_MIN)

# collapse across seeds
agg = defaultdict(list)
for r in rows:
    agg[r["name"]].append(r)

out = []
for name, rs in agg.items():
    d = [x["d_iptm"] for x in rs if x["d_iptm"] is not None]
    j = [x["epitope_jaccard"] for x in rs if x["epitope_jaccard"] is not None]
    out.append(dict(
        name=name, n_seeds=len(rs),
        iptm_mean=round(float(np.mean([x["iptm"] for x in rs])), 4),
        d_iptm_median=round(float(np.median(d)), 4) if d else None,
        d_iptm_worst=round(float(np.min(d)), 4) if d else None,
        win_rate=round(sum(1 for x in d if x > 0) / len(d), 2) if d else None,
        epitope_jaccard_median=round(float(np.median(j)), 3) if j else None,
        pose_ok_rate=round(sum(1 for x in rs if x["pose_ok"]) / len(rs), 2),
        plddt_mean=round(float(np.mean([x["plddt"] for x in rs])), 4),
        ipsae=next((x["ipsae"] for x in rs if x["ipsae"] is not None), None),
    ))

src = {}
if os.path.exists("results/merged60.csv"):
    for r in csv.DictReader(open("results/merged60.csv")):
        src[r["name"]] = (r["source"], r["n_mut"], r["n_zoneA"], r["mutations"])
for o in out:
    s = src.get(o["name"], ("control", "", "", ""))
    o["source"], o["n_mut"], o["n_zoneA"], o["mutations"] = s

out.sort(key=lambda x: (x["d_iptm_median"] is None, -(x["d_iptm_median"] or 0)))
cols = ["name", "source", "n_seeds", "d_iptm_median", "d_iptm_worst", "win_rate",
        "iptm_mean", "epitope_jaccard_median", "pose_ok_rate", "plddt_mean",
        "ipsae", "n_mut", "n_zoneA", "mutations"]
with open("results/screen_scores.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
    w.writeheader()
    w.writerows(out)

print(f"scored {len(rows)} predictions over {len(agg)} sequences\n")
dec = [o for o in out if o["name"].startswith("decoy")]
par = [o for o in out if o["name"] == "parent"]
print("CONTROLS")
for o in par + dec:
    print(f"  {o['name']:20} iptm {o['iptm_mean']:.4f}  d {str(o['d_iptm_median']):>8}  "
          f"jaccard {o['epitope_jaccard_median']}  pose_ok {o['pose_ok_rate']}")
band = max((abs(o["d_iptm_median"]) for o in dec if o["d_iptm_median"] is not None), default=0.04)
print(f"\nnoise band from decoys: +/- {band:.4f} -- treat |d_iptm| below this as indistinguishable\n")

print("TOP 15 by median paired delta")
print(f"{'name':26} {'src':8} {'d_iptm':>8} {'worst':>8} {'win':>5} {'jac':>6} {'pose':>5}")
for o in out[:15]:
    if o["name"] in ("parent",) or o["name"].startswith("decoy"):
        continue
    print(f"{o['name']:26} {o['source']:8} {str(o['d_iptm_median']):>8} "
          f"{str(o['d_iptm_worst']):>8} {str(o['win_rate']):>5} "
          f"{str(o['epitope_jaccard_median']):>6} {str(o['pose_ok_rate']):>5}")
print("\nwrote results/screen_scores.csv")
