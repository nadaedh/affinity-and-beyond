# HER2 final-15 v1

This is an evidence-balanced laboratory portfolio, not a claim of experimentally improved challenge affinity.

## Portfolio composition

- 7 computational leaders: positive median paired ipSAE change, with 75–100% parent-like pose recovery across four seeds.
- 4 structural backups: chosen for pose stability, matched-pair value, and mutation diversity rather than raw ipSAE alone.
- 4 experimental anchors: trastuzumab-scaffold HCDR3 variants with independent published SPR measurements. Their KD values are not directly comparable with the challenge parent because assay/platform conditions may differ.

## Current leader

`LL_CORE_TTK` is the strongest computational result. It won against the paired parent in all four seeds and was the only main-screen candidate with a positive worst-case paired ipSAE change.

## Focused round-two result

`R2_A72R_TTK` passed four-seed confirmation: median paired delta +0.024274, worst-case delta -0.034543, 75% wins, 100% parent-like pose recovery, and median epitope Jaccard 0.880.

## Exclusions

Known unstable or misdocked designs such as `RAT_02`, the L55E branch, and the failing T93R/A72R backgrounds were not used to pad the portfolio.

## Validation

The included FASTA has exactly 15 unique sequences. Every design has at least five edits, is at most 260 aa, contains one intact `(GGGGS)3` linker, retains two VH and two VL conserved cysteines without adding cysteines, and uses canonical amino acids only.

The exact organiser submission schema should still be checked before upload. Computational scores do not establish improved affinity; laboratory measurements are required.
