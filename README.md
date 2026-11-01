# Affinity and Beyond

**Optimising the trastuzumab scFv against HER2**

Breast Cancer Awareness Hackathon, Zurich, 31 October 2026. Protein design
track, run by [Adaptyv Bio](https://adaptyvbio.com) and HerCode.

---

## The challenge

Take trastuzumab (a breast cancer antibody), as a 242 aa single-chain
variable fragment, and make it bind HER2 more tightly.

| | |
|---|---|
| target | HER2 extracellular domain, UniProt P04626 |
| structure | [PDB 1N8Z](https://www.rcsb.org/structure/1N8Z), HER2 + trastuzumab Fab, 2.52 Å |
| parent | 242 aa, V<sub>H</sub>–(GGGGS)<sub>3</sub>–V<sub>L</sub>, **K<sub>D</sub> = 0.82 nM** by BLI |
| rules | ≥5 mutations, one chain, ≤260 aa |
| prize | 15 designs made and measured in Adaptyv's automated lab |

The starting point is not a weak binder. It is an approved drug that has
already been through industrial affinity maturation.

## What we submitted

**[`submission/final15.csv`](submission/final15.csv)** — 15 designs, each 5–6
mutations from the parent. All pass the rules; checked independently in
[`submission/VALIDATION.md`](submission/VALIDATION.md).

They are being expressed cell-free and measured by BLI in Lausanne. Results
land around 21 November 2026, public on [ProteinBase](https://proteinbase.com),
and will be recorded in [`RESULTS.md`](RESULTS.md).

## The result worth reading

We screened **160 variants with 10 samples each** and scored them with ipSAE,
the interface confidence metric the challenge is judged on.

**Not one design is distinguishable from the parent.**

| | |
|---|---|
| parent ipSAE_min | 0.7854 ± 0.0372 |
| best variant (TB_030) | 0.8209 ± 0.0371 |
| difference | +0.0355, t = 2.14, p = 0.033 |
| Bonferroni threshold for 160 tests | p < 0.000313 |
| **designs significant after correction** | **0 / 160** |

The best p-value in the entire screen is **58× too large** to survive
correction for the number of comparisons made.

### The metric disagrees with itself

The **same parent sequence**, scored by two independent Boltz pipelines built
by two people on the same day:

```
Dominik  (10 samples)  0.7854
nada     ( 1 seed)     0.6448
                       ------
gap                    0.1406   <- 4.3x the best variant effect
```

A measurement that varies more between implementations than between molecules
cannot rank molecules.

Two further checks agree:

- **Rank correlation across pipelines: r = −0.322.** Ten of Saniya's designs
  were scored by both her pipeline and nada's. The two runs disagree about
  which are good — the correlation is *negative*.
- **Sign flips within one pipeline: 6 of 11.** More than half of Saniya's
  four-seed designs change sign between their median and their worst seed.

### We calibrated with deliberately broken controls

Rather than trusting the score, we scored designs we knew were damaged:

| control | ipSAE_min | vs parent |
|---|---|---|
| parent | 0.6448 | — |
| Y105A — delete the largest anchor (79 Å² buried) | 0.5570 | −0.088 |
| scrambled CDR-H3 — destroy the loop | 0.4153 | −0.230 |

Useful: ipSAE separates *broken* from *intact*. Not useful: it cannot separate
*good* from *slightly better*. So we used it as a filter, not a ranking —
which is what the organisers advised.

## How we got there

Five people worked in parallel with different methods, and that turned out to
be the most interesting part.

1. **Map the structure first.** From 1N8Z we computed which residues actually
   touch HER2: 20 direct contacts under 4.5 Å, 47 second-shell, the rest
   framework. As an independent check, **all 20 contact residues fall inside
   CDRs** defined purely by sequence.

2. **Point each model only where it is competent.** AbLang2 (antibody
   repertoires), ESM2 (general proteins) and AbMPNN (antibody structures) were
   made to veto each other. This caught real errors — ESM2 scored a
   solvent-exposed framework phenylalanine at **+9.9** where AbLang2 scored it
   **−9.2**. The antibody model was right.

3. **Both design models wanted to destroy the antibody.** AbLang2 and AbMPNN
   share no training data, and both independently proposed reverting the
   paratope. AbMPNN keeps the parent residue **0% of the time** at the three
   largest anchors. They are correct that trastuzumab's binding site is
   atypical — and wrong about what to do, because that atypicality is what
   affinity maturation *produces*. Typicality is what got optimised away.

4. **Fix liabilities structurally, not cosmetically.** The parent carries an
   Asn-Gly deamidation motif and an Asp-Gly isomerisation motif. The
   deamidation site sits in the second shell and is free to repair. The
   isomerisation site is adjacent to **G103, which buries 63 Å² across four
   HER2 contacts** — so we repaired at the aspartate instead and left the
   anchor intact. Three independent pools reached the same conclusion.

### Four methods, one answer

No coordination, four approaches, and the same mutations kept appearing:

| mutation | found independently by |
|---|---|
| **A72R** | nada, saniya, annina-max — **18 designs** |
| N55S / N55Q | nada, saniya, annina-max |
| M107L / M107F | nada, dominik, annina-max |
| D102E / D102S | nada, dominik, annina-max |
| S97A, V168L, A79V, T54D | nada, annina-max |

Across all pools the team explored **42 distinct positions** of 227. Two of us
independently built a **scrambled-CDR-H3 negative control** without discussing
it.

## What we would do differently

- **Replicates before rankings.** A single seed cannot resolve differences
  this small. Dominik's 10-sample screen is the only part of this project that
  supports a statistical claim, and the claim it supports is a null one.
- **Diversify harder.** The submitted 15 carries T32Y in 11 designs and M107L
  in 10, spanning 20 positions. When the ranking is not statistically
  meaningful, the portfolio should be *more* spread out, not less. An
  alternative 15 covering 32 positions sits unused in
  [`methods/nada/nada_pool_final15.fasta`](methods/nada/nada_pool_final15.fasta).
- **Treat the metric as the experiment.** The most informative thing we did was
  score things we knew were broken.

## Repository layout

```
submission/     the 15 designs sent to Adaptyv, and their validation
analysis/       the screens, and the script that reproduces every number here
methods/        each contributor's pipeline, kept separate on purpose
data/           1N8Z, the parent scFv, the trimmed HER2 domain IV target
RESULTS.md      the measured KD values, when they arrive
```

Every figure quoted above is regenerated by:

```bash
python analysis/reproduce_analysis.py
```

It needs only the Python standard library and the CSVs in `analysis/screens/`.

## The open question

In three weeks Adaptyv returns K<sub>D</sub>, k<sub>on</sub>, k<sub>off</sub>
and raw sensorgrams for all 15 designs.

We have published, in advance, the prediction that **our computational scores
carry no information about which of them binds better**. If a design beats
0.82 nM, we will not be able to claim we picked it for the right reason. If
none does, the null result above will have been honest.

Either way the answer goes in [`RESULTS.md`](RESULTS.md).

## Contributors

Nada Eddahir · Sanny · Dominik · Annina · Stella

Built in one day using [Adaptyv's protein design skills](https://github.com/adaptyvbio/protein-design-skills)
for Claude Code, Boltz-2 on Modal GPUs, AbLang2, ESM2 and AbMPNN.

## Licence

Code MIT, data CC BY 4.0 — see [LICENSE](LICENSE). PDB 1N8Z is from the RCSB
and carries its own terms. AbMPNN weights are CC BY 4.0 (Exscientia).
