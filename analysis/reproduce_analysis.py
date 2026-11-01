"""Reproduce every number quoted in the analysis write-ups.

    python analysis/reproduce_analysis.py

Needs only numpy-free stdlib plus the CSVs in analysis/screens/.
"""
import csv
import itertools
import math
import os
import statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
PARENT = ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
          "GGGGSGGGGSGGGGS"
          "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK")


def rd(p):
    with open(os.path.join(ROOT, p)) as fh:
        return list(csv.DictReader(fh))


def norm(p):
    """two-sided p from a z/t statistic"""
    return 2 * (1 - 0.5 * (1 + math.erf(abs(p) / math.sqrt(2))))


def hdr(t):
    print("\n" + "=" * 68)
    print(t)
    print("=" * 68)


# ---------------------------------------------------------------- screen 1
hdr("1. Dominik's library screen: is ANY variant above the parent?")
rows = rd("analysis/screens/dominik_library_screen_ipsae.csv")
p = next(r for r in rows if r["name"].startswith("PARENT"))
pm, ps, pn = float(p["ipSAE_min"]), float(p["ipSAE_sd"]), int(p["n_samples"])
var = [r for r in rows if not r["name"].startswith("PARENT")]
print(f"variants: {len(var)}   parent ipSAE_min = {pm:.4f} +/- {ps:.4f} (n={pn})")

res = []
for r in var:
    m, s, n = float(r["ipSAE_min"]), float(r["ipSAE_sd"]), int(r["n_samples"])
    se = math.sqrt(ps ** 2 / pn + s ** 2 / n)
    t = (m - pm) / se if se else 0.0
    res.append((t, m - pm, r["name"], m, s))
res.sort(reverse=True)

above = sum(1 for r in var if float(r["ipSAE_min"]) > pm)
print(f"above the parent mean: {above}/{len(var)} ({above/len(var):.0%})")
print(f"\n{'variant':10} {'ipSAE':>7} {'sd':>7} {'delta':>8} {'t':>6} {'p':>9}")
for t, d, n, m, s in res[:5]:
    print(f"{n:10} {m:7.4f} {s:7.4f} {d:+8.4f} {t:6.2f} {norm(t):9.5f}")

bonf = 0.05 / len(var)
sig = [x for x in res if x[1] > 0 and norm(x[0]) < bonf]
print(f"\nBonferroni threshold for {len(var)} tests: p < {bonf:.6f}")
print(f"variants significantly above parent: {len(sig)}")
print(f"best p in the whole screen: {norm(res[0][0]):.5f}  "
      f"({norm(res[0][0])/bonf:.0f}x above the threshold)")

# ---------------------------------------------------------------- screen 2
hdr("2. The same parent, scored by two pipelines")
nada = rd("analysis/screens/nada_merged60_screen.csv")
npar = next((r for r in nada if r["name"] == "parent"), None)
if npar and npar.get("ipsae"):
    nv = float(npar["ipsae"])
    print(f"Dominik  (Boltz, 10 samples) : {pm:.4f}")
    print(f"nada     (Boltz,  1 seed)    : {nv:.4f}")
    gap = abs(pm - nv)
    best = res[0][1]
    print(f"\ngap between pipelines        : {gap:.4f}")
    print(f"best variant effect (Dominik): {best:.4f}")
    print(f"ratio                        : {gap/best:.1f}x")
    print("\nThe same molecule scores further apart on two implementations")
    print("than the best design scores from the parent.")

# decoys from nada's screen
dec = [r for r in nada if r["name"].startswith("decoy")]
if dec and npar:
    print(f"\ndecoy calibration (nada pipeline, ipSAE_min):")
    print(f"  parent            {float(npar['ipsae']):.4f}")
    for d in dec:
        if d.get("ipsae"):
            print(f"  {d['name']:17} {float(d['ipsae']):.4f}  "
                  f"({float(d['ipsae'])-float(npar['ipsae']):+.4f})")

# ---------------------------------------------------------------- screen 3
hdr("3. Cross-pipeline rank agreement on shared designs")
san = rd("analysis/screens/saniya_4seed_selection.csv")
hers = {r["candidate_id"]: float(r["paired_delta_median"])
        for r in san if r.get("paired_delta_median")}
ours = {r["name"][4:]: float(r["d_iptm_median"]) for r in nada
        if r["name"].startswith("san_") and r.get("d_iptm_median")}
pairs = [(hers[k], ours[k]) for k in hers if k in ours]
if len(pairs) > 2:
    xs, ys = [a for a, _ in pairs], [b for _, b in pairs]
    mx, my = st.mean(xs), st.mean(ys)
    num = sum((a - mx) * (b - my) for a, b in pairs)
    den = math.sqrt(sum((a - mx) ** 2 for a in xs) * sum((b - my) ** 2 for b in ys))
    print(f"designs scored by both pipelines: {len(pairs)}")
    print(f"Pearson r between her 4-seed medians and our 1-seed deltas: {num/den:+.3f}")
flips = sum(1 for r in san
            if r.get("paired_delta_median") and r.get("paired_delta_worst")
            and float(r["paired_delta_median"]) > 0 > float(r["paired_delta_worst"]))
tot = sum(1 for r in san if r.get("paired_delta_median"))
print(f"designs whose sign flips between median and worst seed: {flips}/{tot}")

# ---------------------------------------------------------------- submission
hdr("4. The submitted 15: legality and diversity")
sub = rd("submission/final15.csv")
bad = 0
for r in sub:
    s = r["sequence"]
    edits = sum(1 for a, b in zip(PARENT, s) if a != b)
    ok = (len(s) <= 260 and edits >= 5 and s.count("GGGGSGGGGSGGGGS") == 1
          and s.count("C") == 4 and set(s) <= set("ACDEFGHIKLMNPQRSTVWY"))
    if not ok:
        bad += 1
        print(f"  FAIL {r['candidate_id']}")
print(f"{len(sub)} designs, {len({r['sequence'] for r in sub})} unique, "
      f"{bad} rule violations")

freq = {}
for r in sub:
    for m in r["edits_1_based_scFv"].split(";"):
        freq[m.strip()] = freq.get(m.strip(), 0) + 1
print("\nmost common mutations:")
for k, c in sorted(freq.items(), key=lambda x: -x[1])[:6]:
    print(f"  {k:8} {c:2}/{len(sub)}")
pos = sorted({int(k[1:-1]) for k in freq})
print(f"\ndistinct positions: {len(pos)}")
sets = [set(r["edits_1_based_scFv"].split(";")) for r in sub]
j = [len(a & b) / len(a | b) for a, b in itertools.combinations(sets, 2)]
print(f"mean pairwise Jaccard: {st.mean(j):.3f}   max: {max(j):.3f}")
print(f"pairs sharing >=3 of 5 mutations: "
      f"{sum(1 for a, b in itertools.combinations(sets, 2) if len(a & b) >= 3)}")

hdr("done")
