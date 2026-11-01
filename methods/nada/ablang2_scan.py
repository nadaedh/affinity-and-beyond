"""AbLang2 paired masked scan over the trastuzumab scFv.

Scores every single substitution as  log P(mut) - log P(wt)  at that position,
from the paired VH/VL model, and annotates each with its interface zone.

Writes results/ablang2_scan.csv
"""
import csv
import os
import sys

import numpy as np
import torch
import ablang2

VH = "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
VL = "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
VH_LEN, LINKER_LEN = 120, 15
AAS = set("ACDEFGHIKLMNPQRSTVWY")

# zones, as computed by scripts/zones_and_liabilities.py (scFv numbering)
ZONE_A = {33, 50, 52, 57, 58, 59, 99, 101, 102, 103, 104, 105,
          163, 165, 166, 167, 184, 185, 188, 201, 226, 227, 228, 229}
ZONE_B = {29, 30, 31, 32, 34, 35, 47, 48, 49, 51, 53, 54, 55, 56, 60, 61, 70,
          97, 98, 100, 106, 107, 108,
          137, 161, 162, 164, 168, 169, 181, 182, 183, 186, 187, 189, 190,
          203, 204, 206, 224, 225, 230, 231}
# never propose changes here: conserved disulfide cysteines
FORBIDDEN = {22, 96, 158, 223}


def zone(i):
    return "A" if i in ZONE_A else ("B" if i in ZONE_B else "C")


def scfv_pos(domain, local):
    """local is 1-based within VH or VL -> position in the 242 aa scFv"""
    return local if domain == "VH" else local + VH_LEN + LINKER_LEN


model = ablang2.pretrained("ablang2-paired", device="cpu")
tok = model.tokenizer

logits = np.asarray(model([[VH, VL]], mode="likelihood")[0])

# --- establish token -> residue mapping, then VERIFY it ---
encoded = tok([[VH, VL]], pad=True, w_extra_tkns=True, device="cpu")
ids = [int(i) for i in encoded[0]]
inv = {v: k for k, v in tok.aa_to_token.items()}
toks = [inv.get(i, "?") for i in ids]

residue_idx = [n for n, t in enumerate(toks) if t in AAS]
if len(residue_idx) != len(VH) + len(VL):
    sys.exit(f"token mapping failed: {len(residue_idx)} residues vs {len(VH)+len(VL)} expected")

chain_seq = VH + VL
for n, ti in enumerate(residue_idx):
    if toks[ti] != chain_seq[n]:
        sys.exit(f"token mismatch at {n}: {toks[ti]} vs {chain_seq[n]}")
print(f"token mapping verified: {len(residue_idx)} residues, logits {logits.shape}")

logprobs = torch.log_softmax(torch.tensor(logits), dim=-1).numpy()

rows = []
for n, ti in enumerate(residue_idx):
    if n < len(VH):
        domain, local = "VH", n + 1
    else:
        domain, local = "VL", n - len(VH) + 1
    i = scfv_pos(domain, local)
    wt = chain_seq[n]
    if i in FORBIDDEN:
        continue
    wt_lp = logprobs[ti, tok.aa_to_token[wt]]
    for mut in sorted(AAS):
        if mut == wt:
            continue
        score = float(logprobs[ti, tok.aa_to_token[mut]] - wt_lp)
        rows.append(dict(scfv_pos=i, domain=domain, local=local, wt=wt,
                         mut=mut, zone=zone(i), ablang2=round(score, 4)))

rows.sort(key=lambda r: -r["ablang2"])
os.makedirs("results", exist_ok=True)
with open("results/ablang2_scan.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)

print(f"wrote results/ablang2_scan.csv  ({len(rows)} substitutions)")
print(f"\npositive-scoring substitutions: {sum(r['ablang2'] > 0 for r in rows)}")
for z in "ABC":
    sub = [r for r in rows if r["zone"] == z and r["ablang2"] > 0]
    print(f"  zone {z}: {len(sub)}")

print("\ntop 25 overall:")
print(f"{'pos':>5} {'dom':4} {'mut':8} {'zone':5} {'score':>7}")
for r in rows[:25]:
    print(f"{r['scfv_pos']:5} {r['domain']:4} {r['wt']}->{r['mut']:<5} {r['zone']:5} {r['ablang2']:7.3f}")
