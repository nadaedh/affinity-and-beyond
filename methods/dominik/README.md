# dominik

A 160-variant combinatorial library, each scored with **10 Boltz samples**.

| file | what |
|---|---|
| `sequences_ipsae.csv` | 160 variants + parent: ipSAE_min, **standard deviation over 10 samples**, epitope recall, mutation list and full sequence |

## Why this is the most important file in the repository

It is the only screen here with replicates, so it is the only one that
supports a statistical test rather than a ranking. That test is what produces
the project's main finding:

**0 of 160 designs are significantly above the parent** after correcting for
160 comparisons. Best variant t = 2.14, p = 0.033, Bonferroni threshold
p < 0.000313.

See [`../../analysis/README.md`](../../analysis/README.md) section 1.

13 of the 15 submitted designs come from this screen.

## Missing

The scripts that generated it were not committed to any branch, so the screen
cannot be re-run from this repository — only re-analysed. The analysis is in
[`../../analysis/reproduce_analysis.py`](../../analysis/reproduce_analysis.py).
