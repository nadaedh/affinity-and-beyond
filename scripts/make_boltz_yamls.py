"""Write one Boltz YAML per sequence: chain A = scFv, chain B = HER2 domain IV.

Outputs results/boltz_inputs/<name>.yaml
"""
import os

HER2 = open("data/her2_domainIV.fasta").read().split("\n")[1].strip()
OUT = "results/boltz_inputs"
os.makedirs(OUT, exist_ok=True)

seqs, name = {}, None
for line in open("results/merged60.fasta"):
    line = line.strip()
    if line.startswith(">"):
        name = line[1:].split()[0]
        seqs[name] = ""
    elif line:
        seqs[name] += line

TEMPLATE = """sequences:
    - protein:
        id: A
        sequence: {scfv}
    - protein:
        id: B
        sequence: {her2}
"""

for n, s in seqs.items():
    with open(f"{OUT}/{n}.yaml", "w") as fh:
        fh.write(TEMPLATE.format(scfv=s, her2=HER2))

print(f"wrote {len(seqs)} YAMLs to {OUT}/")
print(f"scFv length {len(next(iter(seqs.values())))}, HER2 domain IV {len(HER2)}")
print(f"complex size {len(next(iter(seqs.values()))) + len(HER2)} aa")
print("\nfirst 5:", sorted(seqs)[:5])
