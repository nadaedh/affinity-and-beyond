# Measured results

**Status: awaiting lab data.**

The 15 designs in [`submission/final15.csv`](submission/final15.csv) are being
made and measured by Adaptyv Bio in Lausanne: cell-free expression, then
binding kinetics by BLI against the HER2 extracellular domain.

Expected around **21 November 2026**, published on
[ProteinBase](https://proteinbase.com).

## The prediction we are on record for

Before seeing any measurement, this project concluded that its own
computational scores **cannot rank these designs**. See the README and
[`analysis/`](analysis/) for the evidence:

- 0 of 160 variants significantly above parent after Bonferroni correction
- the same parent sequence scores 0.7854 and 0.6448 on two pipelines, a gap
  4.3x larger than the best variant effect
- cross-pipeline rank correlation r = -0.322 on shared designs

So the honest prior is: **no relationship between our ipSAE ranking and the
measured KD.**

## Awaiting results:

| | |
|---|---|
| designs expressed | / 15 |
| designs with measurable binding | / 15 |
| best KD | nM |
| parent KD (reference) | 0.82 nM |
| designs beating parent | |
| Spearman r, ipSAE vs measured KD | |

## If a design does beat 0.82 nM

It should be said plainly: we would not be able to claim we selected it for
the right reason. The selection scores do not separate it from the others.
