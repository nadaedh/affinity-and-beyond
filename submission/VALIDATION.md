# Validation of the submitted 15

Checked independently of the pipeline that produced the designs, by
`analysis/reproduce_analysis.py` section 4.

## Challenge rules

| rule | result |
|---|---|
| ≥ 5 mutations from the parent | pass, all 15 (13 have 5, one has 6) |
| one chain, V<sub>H</sub>–(GGGGS)<sub>3</sub>–V<sub>L</sub> | pass, exactly one intact linker each |
| ≤ 260 aa | pass, all 242 aa |
| canonical amino acids only | pass |
| conserved cysteines retained, none added | pass, exactly 4 each (C22, C96, C158, C223) |
| unique sequences | pass, 15 / 15 |

**0 rule violations.**

## Provenance

| slots | source | selection basis |
|---|---|---|
| 1–11, 14, 15 | Dominik's 160-variant screen | ipSAE_min over 10 samples, epitope recall |
| 12–13 | Saniya's campaign | paired 4-seed delta against parent |

Verified: each design's stated mutation list matches its sequence when
compared position-by-position against the parent.

## Caveats stated up front

1. The file as committed carries `submission_status =
   provisional_needs_rules_and_final_review` for every row. It is recorded
   here as submitted; if the final upload differed, this file should be
   corrected.

2. **The selection scores are not statistically significant.** None of the 160
   screened variants separates from the parent after correcting for the number
   of comparisons. See [`../analysis/README.md`](../analysis/README.md)
   sections 1–3. The 15 are a defensible, rule-compliant, diverse-ish
   portfolio — they are not a ranked list of likely improvements, and should
   not be read as one.

3. Mutation positions are 1-based in the 242 aa scFv (V<sub>H</sub> 1–120,
   linker 121–135, V<sub>L</sub> 136–242). **Not Kabat.** Literature positions
   must be mapped before comparison.
