"""Compare the three team variant pools: ours, Saniya's final 15, the team library."""
import csv
import sys
from collections import Counter

PARENT = ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
          "GGGGSGGGGSGGGGS"
          "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK")
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
pools["saniya"] = load_fasta(sys.argv[1])
lib = {}
for r in csv.DictReader(open(sys.argv[2])):
    lib[r["name"]] = r["sequence"]
lib.pop("parent_trastuzumab_scFv", None)
pools["library"] = lib

print("pool sizes:", {k: len(v) for k, v in pools.items()})
print()

# exact sequence overlap
names = list(pools)
for i in range(len(names)):
    for j in range(i + 1, len(names)):
        a, b = names[i], names[j]
        sa, sb = set(pools[a].values()), set(pools[b].values())
        shared = sa & sb
        print(f"exact sequence overlap {a} vs {b}: {len(shared)}")
        for s in list(shared)[:5]:
            na = [k for k, v in pools[a].items() if v == s][0]
            nb = [k for k, v in pools[b].items() if v == s][0]
            print(f"    {na}  ==  {nb}")
print()

mp = {k: {n: muts(s) for n, s in v.items() if muts(s) is not None} for k, v in pools.items()}
for k in mp:
    allpos = sorted({int(m[1:-1]) for s in mp[k].values() for m in s})
    print(f"{k:8} {len(mp[k])} same-length variants, {len(allpos)} distinct positions mutated")
    byzone = Counter(pos[p]["zone"] for p in allpos)
    print(f"         positions: {allpos}")
    print(f"         zones: {dict(byzone)}")
print()

for i in range(len(names)):
    for j in range(i + 1, len(names)):
        a, b = names[i], names[j]
        pa = {int(m[1:-1]) for s in mp[a].values() for m in s}
        pb = {int(m[1:-1]) for s in mp[b].values() for m in s}
        ma = {m for s in mp[a].values() for m in s}
        mb = {m for s in mp[b].values() for m in s}
        print(f"{a} vs {b}: shared positions {sorted(pa & pb)}")
        print(f"{' ':{len(a)+4}}shared exact mutations {sorted(ma & mb) or 'none'}")
print()

# library internal diversity
libmuts = [mp["library"][n] for n in mp["library"]]
sizes = Counter(len(m) for m in libmuts)
print("library edit-distance distribution:", dict(sorted(sizes.items())))
under5 = [n for n, m in mp["library"].items() if len(m) < 5]
print(f"library variants with <5 mutations: {len(under5)}", under5[:10])
dupes = [s for s, c in Counter(pools["library"].values()).items() if c > 1]
print(f"duplicate sequences inside library: {len(dupes)}")
cys = [n for n, s in pools["library"].items() if s.count("C") != 4]
print(f"library variants without exactly 4 cysteines: {len(cys)}", cys[:5])
lnk = [n for n, s in pools["library"].items() if s.count("GGGGSGGGGSGGGGS") != 1]
print(f"library variants without exactly one intact linker: {len(lnk)}", lnk[:5])
