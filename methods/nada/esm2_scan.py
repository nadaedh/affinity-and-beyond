"""ESM2 zero-shot scan over the 242 aa trastuzumab scFv.

Masked-marginal scoring: mask each position in turn, read the model's
distribution there, score every substitution as log P(mut) - log P(wt).
More faithful than the single-pass wt-marginal, and affordable here because
the construct is short.

Writes results/esm2_scan.csv
"""
import csv
import os

import torch
import esm

SCFV = ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
        "GGGGSGGGGSGGGGS"
        "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK")
VH_LEN, LINKER_LEN = 120, 15
LINKER = range(VH_LEN + 1, VH_LEN + LINKER_LEN + 1)   # 121-135, never mutate
AAS = "ACDEFGHIKLMNPQRSTVWY"

ZONE_A = {33, 50, 52, 57, 58, 59, 99, 101, 102, 103, 104, 105,
          163, 165, 166, 167, 184, 185, 188, 201, 226, 227, 228, 229}
ZONE_B = {29, 30, 31, 32, 34, 35, 47, 48, 49, 51, 53, 54, 55, 56, 60, 61, 70,
          97, 98, 100, 106, 107, 108,
          137, 161, 162, 164, 168, 169, 181, 182, 183, 186, 187, 189, 190,
          203, 204, 206, 224, 225, 230, 231}
FORBIDDEN = {22, 96, 158, 223}


def zone(i):
    return "A" if i in ZONE_A else ("B" if i in ZONE_B else "C")


model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
model.eval()
bc = alphabet.get_batch_converter()
_, _, tokens = bc([("scfv", SCFV)])
print(f"model loaded; tokens {tuple(tokens.shape)} for {len(SCFV)} residues")

rows = []
with torch.no_grad():
    for i in range(1, len(SCFV) + 1):
        if i in LINKER or i in FORBIDDEN:
            continue
        masked = tokens.clone()
        masked[0, i] = alphabet.mask_idx          # +0 offset: token 0 is BOS
        logits = model(masked)["logits"]
        lp = torch.log_softmax(logits[0, i], dim=-1)
        wt = SCFV[i - 1]
        wt_lp = lp[alphabet.get_idx(wt)].item()
        for mut in AAS:
            if mut == wt:
                continue
            rows.append(dict(scfv_pos=i, wt=wt, mut=mut, zone=zone(i),
                             esm2=round(lp[alphabet.get_idx(mut)].item() - wt_lp, 4)))
        if i % 40 == 0:
            print(f"  ...{i}/{len(SCFV)}")

rows.sort(key=lambda r: -r["esm2"])
os.makedirs("results", exist_ok=True)
with open("results/esm2_scan.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)

print(f"\nwrote results/esm2_scan.csv ({len(rows)} substitutions)")
pos = [r for r in rows if r["esm2"] > 0]
print(f"positive: {len(pos)}  (A {sum(r['zone']=='A' for r in pos)}, "
      f"B {sum(r['zone']=='B' for r in pos)}, C {sum(r['zone']=='C' for r in pos)})")

print("\ntop 20 in ZONE C (the conservative tier):")
print(f"{'pos':>5} {'mut':9} {'score':>7}")
for r in [r for r in rows if r["zone"] == "C"][:20]:
    print(f"{r['scfv_pos']:5} {r['wt']}->{r['mut']:<6} {r['esm2']:7.3f}")
