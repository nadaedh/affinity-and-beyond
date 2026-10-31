"""Select a merged ~60 variant screening set from all three team pools.

Sources
  ours     results/variant_pool.csv          (53, tier-stratified)
  saniya   final15.fasta on SaniyaSani-patch-1 (11 computational designs)
  library  ~/Downloads/scfv_library_all.csv  (160, framework-heavy)

Picks a quota from each source by greedy max-min diversity over mutation sets,
so the screen covers the team's hypothesis space rather than one pool's.

Outputs results/merged60.csv and results/merged60.fasta
"""
import csv
import os
import sys

PARENT = ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
          "GGGGSGGGGSGGGGS"
          "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK")
QUOTA = {"ours": 22, "saniya": 11, "library": 27}

pos = {int(r["scfv_pos"]): r for r in csv.DictReader(open("results/positions.csv"))}


def muts(s):
    if len(s) != len(PARENT):
        return None
    return frozenset(f"{a}{i+1}{b}" for i, (a, b) in enumerate(zip(PARENT, s)) if a != b)


def load_fasta(p):
    d, n = {}, None
    for line in open(p):
        line = line.strip()
        if line.startswith(">"):
            n = line[1:].split()[0]
            d[n] = ""
        elif line:
            d[n] += line
    return d


pools = {}
pools["ours"] = {r["name"]: r["seq"] for r in csv.DictReader(open("results/variant_pool.csv"))}
pools["saniya"] = {k: v for k, v in load_fasta(sys.argv[1]).items() if not k.startswith("SPR")}
pools["library"] = {r["name"]: r["sequence"]
                    for r in csv.DictReader(open(sys.argv[2]))
                    if not r["name"].startswith("parent")}


def jac(a, b):
    return len(a & b) / len(a | b) if (a | b) else 0.0


def pick(cands, k, already):
    """greedy max-min: repeatedly take the variant least similar to everything chosen"""
    chosen = []
    while len(chosen) < k and cands:
        best, bestscore = None, 2.0
        for n, m in cands.items():
            ref = already + [c[1] for c in chosen]
            worst = max((jac(m, r) for r in ref), default=0.0)
            if worst < bestscore:
                best, bestscore = (n, m), worst
        chosen.append(best)
        cands.pop(best[0])
    return chosen


selected, seen_muts = [], []
for src, k in QUOTA.items():
    cands = {n: muts(s) for n, s in pools[src].items() if muts(s) is not None}
    got = pick(cands, k, seen_muts)
    for n, m in got:
        seen_muts.append(m)
        nA = sum(1 for x in m if pos[int(x[1:-1])]["zone"] == "A")
        selected.append(dict(name=f"{src[:3]}_{n}", source=src, orig=n,
                             n_mut=len(m), n_zoneA=nA,
                             mutations=" ".join(sorted(m, key=lambda x: int(x[1:-1]))),
                             seq=pools[src][n]))
    print(f"{src:8} selected {len(got)}/{k}")

os.makedirs("results", exist_ok=True)
with open("results/merged60.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["name", "source", "orig", "n_mut", "n_zoneA", "mutations", "seq"])
    w.writeheader()
    w.writerows(selected)

# decoys: scrambled CDR-H3, and a known-bad anchor knockout, as negative controls
import random
random.seed(0)
h3 = list(PARENT[98:109])
random.shuffle(h3)
decoy1 = PARENT[:98] + "".join(h3) + PARENT[109:]
decoy2 = PARENT[:104] + "A" + PARENT[105:]       # Y105A, the 79 A^2 anchor removed

with open("results/merged60.fasta", "w") as fh:
    fh.write(f">parent\n{PARENT}\n")
    fh.write(f">decoy_scrambleH3\n{decoy1}\n")
    fh.write(f">decoy_Y105A\n{decoy2}\n")
    for v in selected:
        fh.write(f">{v['name']} {v['mutations']}\n{v['seq']}\n")

print(f"\n{len(selected)} variants + parent + 2 decoys = {len(selected)+3} sequences")
print("wrote results/merged60.csv and results/merged60.fasta")

from collections import Counter
print("\nby source:", dict(Counter(v["source"] for v in selected)))
print("zone-A mutations per variant:", dict(sorted(Counter(v["n_zoneA"] for v in selected).items())))
allpos = sorted({int(m[1:-1]) for v in selected for m in v["mutations"].split()})
print(f"distinct positions covered: {len(allpos)}")
print(allpos)
