# saniya

A focused core-grid campaign with **paired 4-seed deltas** against the parent.

| file | what |
|---|---|
| `final15.fasta` | her 15 designs: 11 computational, 4 published SPR-measured anchors |
| `approach.md` | her own write-up of the portfolio logic |
| `VALIDATION.txt` | her rule checks |

2 of these designs (`LL_CORE_TTK`, `LL_A72R_TYK`) are in the team submission.

## What this approach contributed

The three metrics the rest of the team adopted:

- **paired delta on matched seeds** — compare variant and parent on the *same*
  seed, so seed noise cancels
- **pose stability** — does the predicted pose stay parent-like across seeds
- **epitope Jaccard** — is it still binding the *same site*, which catches a
  design that scores well while docking elsewhere

## What her own data shows

**6 of 11 designs flip sign** between their median and their worst seed. Her
`LL_CORE_SYK` has a median of +0.095 and a worst case of −0.123.

That is the earliest evidence in this project that single-run rankings are
unreliable — visible inside one pipeline, before any cross-pipeline
comparison was made.

The 4 SPR anchors are published trastuzumab-scaffold variants with measured
KD (0.63, 0.64, 0.87, 0.94 nM). Two are **weaker** than the 0.82 nM parent,
and the platform differs, so they are reference points rather than candidates.
