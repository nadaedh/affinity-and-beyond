# Methods

Four approaches, run in parallel by different people on the same day, with no
shared pipeline. They are kept separate on purpose: the agreement between them
is evidence precisely *because* they are independent. Merging them into one
tidy pipeline would destroy that.

| directory | approach | output |
|---|---|---|
| [`nada/`](nada/) | structure-derived interface zones; AbLang2 x ESM2 x AbMPNN with mutual veto; Boltz-2 screen of 60 merged variants with two decoy controls | 15 designs spanning 32 positions (not submitted) |
| [`dominik/`](dominik/) | 160-variant combinatorial library, **10 Boltz samples per sequence** with standard deviations and epitope recall | 13 of the 15 submitted designs |
| [`saniya/`](saniya/) | focused core-grid campaign, paired 4-seed deltas against parent, pose and epitope-overlap checks, round-two confirmation | 2 of the 15 submitted designs |
| [`annina-max/`](annina-max/) | interface analysis from 1N8Z, hypothesis-driven variants (safe / disulfide / deamidation / oxidation / charge / paratope), each with its own prediction log, plus a scrambled-H3 control | 15 designs, not submitted |
| [`stella/`](stella/) | orthogonal route: the published Mason et al. deep-learning CDR-H3 library, filtered on charge, GRAVY and pI | 14 designs, not submitted |

## Notes on each

**nada** — the only pool with deliberate negative controls scored alongside
the designs, which is what produced the decoy calibration. Single seed, so
individually underpowered. Full decision log in
[`nada/lab-notebook.md`](nada/lab-notebook.md).

**dominik** — the statistically strongest work here. Ten samples per sequence
is what makes the null result in the README provable rather than suspected.
Only the output table survives; the generating scripts were not committed.

**saniya** — introduced the three metrics the rest of us adopted: paired
deltas on matched seeds, pose stability, and epitope Jaccard against the
parent. Her own data shows 6 of 11 designs flipping sign across seeds.

**annina-max** — the clearest hypothesis-per-design structure: each variant
has a stated mechanism and its own log file. Independently built a scrambled
CDR-H3 control, as nada did, without discussion.

**stella** — the only approach not derived from the 1N8Z interface map. Uses
published, experimentally validated CDR-H3 variants instead. `Top_15_*.csv`
as originally committed had a CSV quoting fault that shifted columns; the
repaired copy is `Top_15_scFv_Varianten_fixed.csv` (14 rows, all 242 aa) and
the original is kept beside it unchanged.

## Reproducing nada's pipeline

```bash
pip install modal ablang2 fair-esm biopython
modal setup
git clone https://github.com/hgbrian/biomodals
git clone --depth 1 https://github.com/dauparas/ProteinMPNN tools/ProteinMPNN
curl -L -o tools/ProteinMPNN/abmpnn_weights/abmpnn.pt \
  https://zenodo.org/records/8164693/files/abmpnn.pt
git clone --depth 1 https://github.com/DunbrackLab/IPSAE.git tools/IPSAE
```

Then, in order: `ablang2_scan.py`, `esm2_scan.py`, `run_abmpnn.sh`,
`parse_abmpnn.py`, `annotate_positions.py`, `build_variants.py`,
`prep_complex.py`, `make_boltz_yamls.py`, `run_boltz_screen.sh`,
`run_ipsae.sh`, `score_screen.py`, `select_final15.py`.

Paths assume the project root; the scripts were written for a working tree
that also contained `biomodals/` and `tools/`.
