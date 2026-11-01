# annina-max

Hypothesis-driven design: each variant states a mechanism and was meant to be
tested on its own.

| file | what |
|---|---|
| `scan_results.csv` | **the valuable one.** 4,314 single mutations scored, with IMGT numbering, CDR/FW assignment, min distance to HER2, shell classification, VH/VL interface flag, AbLang2 and ESM2 scores |
| `1N8Z_fab_interface.csv` | interface residues derived from the crystal structure |
| `her2_scfv_15_variants.csv` / `.fasta` | the 15 designs, grouped by hypothesis (safe / disulfide / deamidation / oxidation / charge / germline / paratope) |
| `prediction_runs.md` | consolidated from 16 Modal logs — **13 of 16 runs failed on rate limiting**, including the control |
| `make_target.py`, `make_boltz_inputs.py`, `build_variants.py`, `score_all.py` | the pipeline |
| `her2_target*.fasta`, `her2_target.pdb`, `parent.json` | target preparation |

`1n8z.pdb` was removed as a duplicate of [`../../data/structures/1n8z.cif`](../../data/structures/1n8z.cif)
(same structure, chains A/B/C, PDB rather than mmCIF format).

## Worth knowing

`scan_results.csv` is an **independent derivation of the same interface map**
nada built, using IMGT numbering and distance shells rather than zones. The
two agree on which positions contact HER2, which is a useful cross-check — see
[`../../analysis/README.md`](../../analysis/README.md).

Because 13 of 16 structure predictions failed, these designs were selected
mainly from the sequence-level scan rather than from cofolding.
