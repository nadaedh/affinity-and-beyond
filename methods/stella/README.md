# stella

The only approach not derived from the 1N8Z interface map. Uses the published
**Mason et al. deep-learning CDR-H3 library** — variants already screened
experimentally — then filters on biophysical properties.

| file | what |
|---|---|
| `mason_validated_scfv.csv` | 30 CDR-H3 variants from the published library, built into the full 242 aa scFv |
| `Top_15_scFv_Varianten_fixed.csv` / `.fasta` | the selection, filtered on charge, GRAVY and pI — **use these** |
| `Top_15_scFv_Varianten.csv` | the original file, kept unchanged (see below) |
| `Top_15_scFv_Varianten_Evaluation.xlsx` | the evaluation spreadsheet |
| `fold.py`, `modal_fold.py`, `scfv_her2.yaml` | Boltz cofolding setup |

## The CSV fault

`Top_15_scFv_Varianten.csv` as originally committed wraps each data row in an
extra pair of quotes with doubled inner quotes around the mutations field. A
standard CSV parser reads the whole row as one field and returns empty
sequences.

The original is kept unchanged; `_fixed.csv` is the repaired copy — 14 rows
(the file is named "Top_15" but holds 14), all 242 aa.

`modal_fold.py` was committed as `import modal.py`; renamed so it can be
imported.

## Why this pool is interesting

Every other approach started from the same crystal structure and the same
interface analysis, so a shared blind spot would be invisible. This one starts
from experimental data instead. None of its designs were submitted, which is
a missed opportunity for diversification given that the selection scores
driving the final 15 are not statistically significant.
