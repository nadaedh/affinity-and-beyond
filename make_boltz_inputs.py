#!/usr/bin/env python3
"""
Write Boltz-2 two-chain input files for the parent, a scrambled control,
and all 15 variants.  Run in ~/her2-scfv after make_target.py.

  python make_boltz_inputs.py

Needs: parent.json, her2_target_full.fasta, her2_scfv_15_variants.fasta
Writes: boltz_in/<name>.fasta   (chain A = scFv, chain B = target)
"""
import json, os, random

TARGET_F = "her2_target_full.fasta"
VARIANTS = "her2_scfv_15_variants.fasta"

def read_fasta(path):
    out, name, buf = {}, None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if name: out[name] = "".join(buf)
            name, buf = line[1:].split("|")[0].strip(), []
        elif line:
            buf.append(line)
    if name: out[name] = "".join(buf)
    return out

target = next(iter(read_fasta(TARGET_F).values()))
assert len(target) == 119, f"target is {len(target)} aa, expected 119"

parent = json.load(open("parent.json"))["full_sequence"]
seqs = {"PARENT": parent}
seqs.update(read_fasta(VARIANTS))

# negative control: scramble CDR-H3 (scFv 97-110), keep length
random.seed(0)
cdr = list(parent[96:110])
while True:
    random.shuffle(cdr)
    if "".join(cdr) != parent[96:110]:
        break
seqs["CTRL_scrambled_H3"] = parent[:96] + "".join(cdr) + parent[110:]

os.makedirs("boltz_in", exist_ok=True)
for name, s in seqs.items():
    assert len(s) == 242, f"{name}: {len(s)} aa"
    assert s[120:135] == "GGGGSGGGGSGGGGS", f"{name}: linker altered"
    with open(f"boltz_in/{name}.fasta", "w") as f:
        f.write(f">A|protein|\n{s}\n>B|protein|\n{target}\n")

print(f"wrote {len(seqs)} input files to boltz_in/")
print(f"  scFv 242 aa + target {len(target)} aa = {242+len(target)} residues per complex")
for n in sorted(seqs):
    print("   ", n)
