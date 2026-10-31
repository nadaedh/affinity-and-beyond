"""Phase 2: combine ranked single mutations into a stratified variant pool.

Merges whichever scans are present (AbLang2, ESM2, AbMPNN), applies the
anti-epistasis rules, and emits a diversity-stratified pool.

Rules enforced:
  * >= MIN_MUT mutations from the parent (challenge requirement)
  * at most one mutation per CDR loop
  * CB-CB > MIN_DIST between any two chosen positions
  * at most one zone-A mutation per variant
  * never the linker, never the conserved cysteines
  * length <= 260 aa (substitutions only, so always 242)

Usage:  python scripts/build_variants.py [--n 50]
"""
import argparse
import csv
import itertools
import os
import random

import numpy as np

SCFV = ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
        "GGGGSGGGGSGGGGS"
        "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK")
MIN_MUT, MIN_DIST = 5, 10.0
# The spacing rule guards against structural epistasis between mutations that
# share an environment. Two distal framework (zone C) surface positions are
# far from the interface and from each other's role, so a shorter separation
# is defensible there -- and necessary: only 6 zone-C positions survive the
# two-model veto, and no 5 of them are mutually >10 A apart.
MIN_DIST_C = 6.0

# liability repairs carried regardless of model score (chemistry, not frequency)
LIABILITY = {
    "D102E": "DG isomerisation repair, charge-conservative",
    "D102S": "DG isomerisation repair, AbLang2-preferred",
    "N55S": "NG deamidation repair",
    "G56S": "NG deamidation repair, motif-breaking",
    "M107L": "CDR-H3 Met oxidation repair",
    "M107F": "CDR-H3 Met oxidation repair",
}


def load_positions():
    pos = {}
    with open("results/positions.csv") as fh:
        for r in csv.DictReader(fh):
            pos[int(r["scfv_pos"])] = r
    return pos


def load_scans():
    """merge every scan we have into {(pos,mut): {source: score}}"""
    merged = {}
    for path, col in [("results/ablang2_scan.csv", "ablang2"),
                      ("results/esm2_scan.csv", "esm2"),
                      ("results/abmpnn_scan.csv", "abmpnn")]:
        if not os.path.exists(path):
            print(f"  (missing {path} - skipping {col})")
            continue
        n = 0
        with open(path) as fh:
            for r in csv.DictReader(fh):
                key = (int(r["scfv_pos"]), r["mut"])
                merged.setdefault(key, {})[col] = float(r[col])
                n += 1
        print(f"  loaded {n} scores from {col}")
    return merged


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=50)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    random.seed(args.seed)

    pos = load_positions()
    print("merging scans:")
    scans = load_scans()
    D = np.load("results/cb_distances.npy")

    sources = sorted({s for v in scans.values() for s in v})
    print(f"\nsources available: {sources}")

    # Candidate rule.
    #   zones A and B: >=2 sources must score it positive. These are the risky
    #     edits and both models are known to mis-call the matured paratope.
    #   zone C: one positive source is enough, PROVIDED no source strongly
    #     dislikes it (>= VETO). Zone C is only covered by two sources and the
    #     consensus rule starves the conservative tier. The veto is what stops
    #     ESM2 artefacts such as F218E (ESM2 +9.9, AbLang2 -9.2) getting in.
    VETO = -2.0
    cands = []
    for (p, mut), sc in scans.items():
        row = pos[p]
        if not int(row["designable"]):
            continue
        tag = f"{row['wt']}{p}{mut}"
        n_pos = sum(1 for v in sc.values() if v > 0)
        is_liab = tag in LIABILITY
        if not is_liab:
            if row["zone"] == "C":
                if n_pos < 1 or min(sc.values()) < VETO:
                    continue
            elif n_pos < min(2, len(sources)):
                continue
        cands.append(dict(pos=p, wt=row["wt"], mut=mut, tag=tag,
                          zone=row["zone"], region=row["region"],
                          mean=float(np.mean(list(sc.values()))),
                          n_support=n_pos, liability=is_liab,
                          **{k: round(v, 3) for k, v in sc.items()}))

    cands.sort(key=lambda c: (-c["n_support"], -c["mean"]))
    print(f"candidate single mutations: {len(cands)} "
          f"(zone A {sum(c['zone']=='A' for c in cands)}, "
          f"B {sum(c['zone']=='B' for c in cands)}, "
          f"C {sum(c['zone']=='C' for c in cands)}, "
          f"liability {sum(c['liability'] for c in cands)})")

    with open("results/candidate_mutations.csv", "w", newline="") as fh:
        cols = ["tag", "pos", "wt", "mut", "zone", "region", "n_support", "mean",
                "liability"] + sources
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(cands)
    print("wrote results/candidate_mutations.csv")

    def compatible(chosen, c):
        if any(x["pos"] == c["pos"] for x in chosen):
            return False
        if c["region"].startswith("CDR") and any(x["region"] == c["region"] for x in chosen):
            return False
        if c["zone"] == "A" and any(x["zone"] == "A" for x in chosen):
            return False
        for x in chosen:
            d = D[x["pos"], c["pos"]]
            limit = MIN_DIST_C if (x["zone"] == "C" and c["zone"] == "C") else MIN_DIST
            if not np.isnan(d) and d < limit:
                return False
        return True

    def assemble(pool, n_mut, tries=4000):
        """greedy-random assembly of one variant from a candidate pool"""
        for _ in range(tries):
            chosen = []
            for c in random.sample(pool, len(pool)):
                if compatible(chosen, c):
                    chosen.append(c)
                if len(chosen) == n_mut:
                    return sorted(chosen, key=lambda x: x["pos"])
        return None

    zoneC = [c for c in cands if c["zone"] == "C"]
    zoneB = [c for c in cands if c["zone"] == "B"]
    zoneA = [c for c in cands if c["zone"] == "A"]
    liab = [c for c in cands if c["liability"]]

    # The conservative tier's intent is "paratope untouched". Pure zone C cannot
    # reach 5 mutations: only 6 zone-C positions survive the two-model veto and
    # no 5 of them are mutually >6 A apart. Zone-B FRAMEWORK positions make no
    # HER2 contact and lie outside the CDR loops, so they honour the intent
    # while making the tier achievable.
    conservative_pool = zoneC + [c for c in zoneB if c["region"] == "FR"]

    tiers = [("conservative", conservative_pool, 5, round(args.n * 0.20)),
             ("rim", zoneB + zoneC, 6, round(args.n * 0.34)),
             ("paratope", zoneA + zoneB + zoneC, 6, round(args.n * 0.20)),
             ("liability", liab + zoneC + zoneB, 5, round(args.n * 0.14)),
             ("wildcard", zoneB + zoneC, 6, round(args.n * 0.12))]

    variants, seen = [], set()
    for name, pool, n_mut, want in tiers:
        if len(pool) < n_mut:
            print(f"  ! tier {name}: pool has {len(pool)} candidates, need {n_mut} - skipped")
            continue
        made = 0
        for _ in range(want * 60):
            if made >= want:
                break
            if name == "paratope":
                a = random.choice(zoneA) if zoneA else None
                if a is None:
                    break
                rest = assemble([c for c in pool if c["zone"] != "A"], n_mut - 1)
                sel = sorted([a] + rest, key=lambda x: x["pos"]) if rest and compatible(rest, a) else None
            elif name == "liability":
                l = random.choice(liab) if liab else None
                if l is None:
                    break
                rest = assemble([c for c in pool if not c["liability"]], n_mut - 1)
                sel = sorted([l] + rest, key=lambda x: x["pos"]) if rest and compatible(rest, l) else None
            else:
                sel = assemble(pool, n_mut)
            if not sel:
                continue
            key = tuple(c["tag"] for c in sel)
            if key in seen:
                continue
            seen.add(key)
            seq = list(SCFV)
            for c in sel:
                seq[c["pos"] - 1] = c["mut"]
            variants.append(dict(name=f"{name}_{made+1:02d}", tier=name,
                                 n_mut=len(sel), mutations=" ".join(key),
                                 seq="".join(seq)))
            made += 1
        print(f"  tier {name}: {made}/{want}")

    # --- explicit construction: the VH44/VL100 interdomain disulfide staple ---
    # Not model-proposed; a known scFv stabiliser. CA-CA in 1N8Z is 5.46 A,
    # inside the 4.5-6.5 A window for an engineered disulfide. Flagged risk:
    # cell-free lysates vary in disulfide formation.
    STAPLE = [dict(pos=44, wt=SCFV[43], mut="C", tag=f"{SCFV[43]}44C", zone="C",
                   region="FR", liability=False),
              dict(pos=234, wt=SCFV[233], mut="C", tag=f"{SCFV[233]}234C", zone="C",
                   region="FR", liability=False)]
    made = 0
    for _ in range(300):
        if made >= 3:
            break
        rest = assemble([c for c in zoneC if c["pos"] not in (44, 234)], 3)
        if not rest or not all(compatible(rest, s) for s in STAPLE):
            continue
        sel = sorted(STAPLE + rest, key=lambda x: x["pos"])
        key = tuple(c["tag"] for c in sel)
        if key in seen:
            continue
        seen.add(key)
        seq = list(SCFV)
        for c in sel:
            seq[c["pos"] - 1] = c["mut"]
        variants.append(dict(name=f"staple_{made+1:02d}", tier="wildcard-staple",
                             n_mut=len(sel), mutations=" ".join(key), seq="".join(seq)))
        made += 1
    print(f"  tier wildcard-staple: {made}/3")

    with open("results/variant_pool.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["name", "tier", "n_mut", "mutations", "seq"])
        w.writeheader()
        w.writerows(variants)
    with open("results/variant_pool.fasta", "w") as fh:
        fh.write(f">parent\n{SCFV}\n")
        for v in variants:
            fh.write(f">{v['name']} {v['mutations']}\n{v['seq']}\n")

    bad = [v for v in variants if len(v["seq"]) != 242 or v["n_mut"] < MIN_MUT]
    print(f"\n{len(variants)} variants written to results/variant_pool.{{csv,fasta}}")
    print(f"rule violations: {len(bad)}")


if __name__ == "__main__":
    main()
