# Lab notebook — trastuzumab scFv optimisation vs HER2

Hackathon: Breast Cancer Awareness, Zurich, 31 Oct 2026.
Target to beat: **0.82 nM KD** (parent VH–(GGGGS)3–VL, BLI, Adaptyv platform).
Deliverable: 15 variants, >=5 mutations each, one chain, <=260 aa.

Numbering convention used throughout: **linear position in the 242 aa scFv**
(VH 1-120, linker 121-135, VL 136-242). *Not* Kabat. Literature that cites
Kabat (e.g. "VH W98") must be mapped before use — see Open questions.

---

## 2026-10-31 — Environment setup

| Step | Result |
|---|---|
| 24 Adaptyv protein-design skills | installed to `.claude/skills/`, gitignored (team already had them) |
| Modal CLI 1.6.3 in `.venv/` | authenticated, profile `nadaedh` |
| `biomodals` cloned | `modal run modal_boltzgen.py --help` verified by hand |

Gotcha: the `setup` skill says to verify auth with `modal token show`.
That subcommand no longer exists in 1.6.3 — use `modal profile current`.

Not yet run: any GPU job. Modal free tier is ~$30/month of credit.

---

## 2026-10-31 — Phase 0 (partial): structure and zones

### Structure
`data/structures/1n8z.cif` from RCSB. Chain content confirmed:

| chain | content | length |
|---|---|---|
| A | VL + CL (light chain) | 214 aa |
| B | VH + CH1 (heavy chain) | 220 aa |
| C | HER2 ectodomain | 581 aa modelled, numbered 1–607 |

VH = 120 aa (**not** 119 — cost us an off-by-one on the first pass).
120 + 15 linker + 107 VL = 242. Parent sequence verified residue-by-residue
against the structure: **0 mismatches**.

### Zone assignment — `scripts/zones_and_liabilities.py`

Of 227 variable-domain residues: **Zone A 20, Zone B 47, Zone C 160.**

**Zone A (direct contact, <4.5 A from HER2)** — 24 residues by the looser
contact-or-burial definition, 20 by strict contact. Ranked by buried area:

| | residues (scFv numbering, dSASA in A^2) |
|---|---|
| VH anchors | Y105 (79), Y57 (65), G103 (63), R59 (62), R50 (38), Y33 (29), W99 (24) |
| VH rim | T58 (11), D102 (8.7), F104 (6.0), Y52 (3.6), G101 (2.4) |
| VL anchors | Y227 (80), N165 (76), F188 (65), T229 (51), H226 (41), T166 (26), S185 (22), R201 (21), D163 (20) |
| VL rim | A167 (12), Y184 (11), T228 (4.7) |

HER2 epitope = residues **557–605**, i.e. domain IV, membrane-proximal.
Consistent with the deck. Supports trimming the target for cofolding.

Worked example of the contact call: VH **Y105** ring stacks on HER2 **F573**
at 3.44 A, and sits at H-bond distance from backbone O of P571 (3.32 A) and
D570 (3.50 A). One residue, three interactions, 79 A^2 buried. Do not touch.

**Zone B (second shell)** — within 5 A of a zone-A side chain, no HER2 contact:
- VH: I29 K30 D31 T32 I34 H35 W47 V48 A49 I51 P53 T54 N55 G56 Y60 A61 I70 S97 R98 G100 G101 D102 A106 M107 D108
- VL: I137 S161 Q162 V164 V168 A169 L181 L182 I183 Y184 A186 S187 L189 Y190 R201 G203 T204 F206 Q224 Q225 P230 P231

### Developability liabilities — checked against zones

| motif | position | zones | risk | verdict |
|---|---|---|---|---|
| **NG** | N55–G56 | B/B | deamidation, fast | **free to fix** — no HER2 contact |
| **DG** | D102–G103 | B/**A** | isomerisation, fast | **partially blocked** — G103 buries 63 A^2 |
| NT | N165–T166 | **A/A** | deamidation, slow | do not touch, N165 buries 76 A^2 |
| NS | N84–S85 | C/C | deamidation | free, low priority |
| NT | N77–T78 | C/C | deamidation, slow | free, low priority |
| DS | D62–S63 | C/C | isomerisation | free, low priority |

**Key call on DG:** the reactive partner in Asp-Gly isomerisation is the
**Asp**, which forms the succinimide. D102 is zone B (8.7 A^2, no direct
contact); G103 is a zone-A anchor. So repair the motif at **D102, not G103**
— D102E removes the succinimide chemistry while leaving the anchor intact.
Mutating G103 to kill the motif would be self-sabotage.

Cysteines: C22, C96, C158, C223 — these are the two canonical intradomain
disulfides (VH 22–96, VL 23–88 local), **not** free thiols. Nothing unpaired.

### Disulfide staple feasibility
Classic scFv VH44/VL100 interdomain staple: in 1N8Z these are Gly (chain B 44)
and Gln (chain A 100), CA–CA **5.46 A**. Engineered disulfides want CA–CA
~4.5–6.5 A, so the geometry works. In scFv numbering that is **G44C + Q234C**,
both zone C.

---

## 2026-10-31 — Phase 1 tooling installed

| tool | how | status |
|---|---|---|
| **AbLang2** | `pip install ablang2` into `.venv` | working on **CPU**. `ablang2.pretrained('ablang2-paired')`, `mode='likelihood'` returns a 232 x 26 logit matrix for the paired VH/VL input (227 residues + separators) |
| **AbMPNN** | weights from [Zenodo 8164693](https://zenodo.org/records/8164693) (`abmpnn.pt`, 20 MB, CC BY 4.0, Exscientia) run with [dauparas/ProteinMPNN](https://github.com/dauparas/ProteinMPNN) code (MIT) | checkpoint verified: `num_edges=48`, `noise_level=0.02` (i.e. the `v_48_020` config), 118 tensors, loads under torch 2.15 |

Both live in `tools/` (gitignored). To reproduce:

```bash
./.venv/bin/pip install ablang2
git clone --depth 1 https://github.com/dauparas/ProteinMPNN.git tools/ProteinMPNN
mkdir -p tools/ProteinMPNN/abmpnn_weights
curl -L -o tools/ProteinMPNN/abmpnn_weights/abmpnn.pt \
  https://zenodo.org/records/8164693/files/abmpnn.pt
```

AbMPNN needs the Fab–HER2 complex as a PDB with the HER2 chain fixed; it is
backbone-fixed design, so it proposes substitutions compatible with the
existing geometry. Runs on CPU for a single complex, slowly.

---

## 2026-10-31 — Organiser guidelines received

Changes to the working plan:

1. **Shortlist of 50–100 is the endorsed cofolding budget**, not 15. Resolves
   the earlier GPU-budget worry — Phase 3 at ~100 complexes is in range.
2. **ipSAE is a filter, not a ranker.** Organisers: "a small ipSAE difference
   between two good variants doesn't mean one binds tighter. Use it to throw
   out broken designs, then pick a diverse 15." So selection = drop broken on
   ipSAE, then **maximise diversity**, not rank-order.
3. **Linker and the four conserved cysteines (C22, C96, C158, C223) are
   off-limits.** Note the VH44/VL100 staple wildcard *adds* a new pair rather
   than touching these — permitted by the letter of the rule, but organisers
   warn cell-free lysates vary in disulfide formation. Cap at 1–2 variants.
4. De novo tools (BindCraft, RFdiffusion, BoltzGen) explicitly discouraged —
   matches the Tier 3 assessment already made. **Boltz-2 ≠ BoltzGen.**
5. Setup should not eat the afternoon. Ours is done.

**Target pool size: ~50, not 15.** Teams submit 15; the track winner gets 45
slots. Build a ranked, diversity-stratified 50 so the winner's expansion set
is ready without a second design round.

---

## 2026-10-31 — Phase 1a: AbLang2 scan (DONE)

`scripts/ablang2_scan.py` -> `results/ablang2_scan.csv`
Score = `log P(mut) - log P(wt)` from the paired model, single forward pass.
Token mapping verified residue-by-residue (227 residues, 232 tokens as
`<VH> | <VL>`). Conserved cysteines excluded from output.

**4237 substitutions scored; only 47 score positive.**
By zone: **A 15, B 28, C 4.**

### Finding 1 — the zone-A hits are germline-reversion traps

AbLang2 is trained on natural repertoires, so "positive" means *more typical
of antibodies in general*, not *better at binding HER2*. Trastuzumab's CDRs
are humanised murine and deliberately atypical. The model's top zone-A calls:

| call | score | what it would destroy |
|---|---|---|
| VL F188S | +1.62 | F188 buries **65 A^2** — one of the four big aromatic anchors |
| VL A167Y | +1.47 | zone-A contact with HER2 571 |
| VL H226S | +1.23 | H226 buries 41 A^2, contacts 571/572 |

Following AbLang2 into zone A would dismantle the paratope. **Confirms the
plan's zone discipline empirically: AbLang2 is for zones B and C only.**

### Finding 2 — M107 oxidation, a third liability the plan missed

VH **M107** sits in CDR-H3 (`WGGDGFYA`**`M`**`DY`), zone B, no HER2 contact.
A methionine in a CDR is a generic oxidation liability. Every hydrophobic
swap is strongly endorsed:

`M107F +2.80`, `M107L +2.12`, `M107I +1.38`, `M107V +1.15`, `M107A +0.77`

Liability repair *and* the model's top structural-neighbourhood call.
(Caveat: documented trastuzumab Met oxidation CQAs are mainly in the **Fc**
(M255/M431), not the Fv — treat the CDR-H3 Met as a plausible liability to
test, not an established one.)

### Finding 3 — zone C is nearly empty, and that is a problem

Only **two** framework positions score positive: VH A72 (R +2.22, K +0.84)
and VH A79 (L +0.42). Expected — trastuzumab is already humanised onto a
human framework, so there is little germline distance left to close.

**Consequence: the conservative tier (zone C only, 5-6 mutations) cannot be
populated from AbLang2 alone.** It needs ESM2 zero-shot or an explicit
germline alignment. This is now the critical gap in Phase 1.

### Finding 4 — least-unnatural liability repairs

| liability | option | AbLang2 | note |
|---|---|---|---|
| NG N55–G56 | N55G | **+0.25** | only positive option; but G-G may over-flexibilise CDR-H2 |
| | G56S | −0.42 | mild; also removes the motif |
| | N55S | −1.29 | standard fix |
| | N55Q / G56A | −4.99 / −3.35 | textbook fixes, strongly unnatural here |
| DG D102–G103 | D102S | **+0.24** | removes the Asp; loses the charge |
| | D102E | −4.47 | charge-conservative, canonical anti-isomerisation fix |
| | G103A | −3.10 | **do not** — G103 is a zone-A anchor (63 A^2) |
| M107 ox | M107F / M107L | **+2.80 / +2.12** | best of both worlds |

**Interpretation note:** a negative AbLang2 score means *rare in natural
repertoires*, not *biophysically bad*. D102E is the canonical isomerisation
repair because Glu cannot form the succinimide; AbLang2 penalises it purely
on rarity. Carry **both D102E and D102S** into the pool and let the lab decide.

### Other notable zone-B calls
`VL Y190Q +3.39` / `Y190E +2.94` (top scorer overall), `VH T32Y +1.31`,
`VH K30S +0.93`, `VH A106Y +0.72`, `VH T54D +0.62`.

---

## 2026-10-31 — Trimmed complex built

`scripts/prep_complex.py`:

| output | content |
|---|---|
| `data/structures/fv_her2_domainIV.pdb` | VH 1-120 (B) + VL 1-107 (A) + HER2 490-607 (C) |
| `data/her2_domainIV.fasta` | **108 aa** HER2 domain IV target for Boltz |
| `results/abmpnn_design_positions.txt` | 35 VH + 32 VL design positions (zones A+B, cysteines excluded) |

HER2 chain C is modelled only to residue 607, so the trim is 490-607, not
490-620. Cofolding complex is **242 + 108 = 350 aa** vs 872 untrimmed.

AbMPNN run script ready at `scripts/run_abmpnn.sh` (user runs it; cloned-code
permission block applies to Claude).

---

## 2026-10-31 — Phase 2 scaffolding: CDR map + distance matrix

`scripts/annotate_positions.py` -> `results/positions.csv`,
`results/cb_distances.npy` (227x227 CB-CB, CA substituted for glycine).

CDRs located by **exact substring match** against the canonical trastuzumab
loops, so no numbering scheme is assumed:

| CDR | sequence | scFv positions | zone A | zone B |
|---|---|---|---|---|
| CDR-H1 | DTYIH | 31–35 | 33 | 31 32 34 35 |
| CDR-H2 | RIYPTNGYTRYADSVKG | 50–66 | 50 52 57 58 59 | 51 53 54 55 56 60 61 |
| CDR-H3 | WGGDGFYAMDY | 99–109 | 99 101 102 103 104 105 | 100 106 107 108 |
| CDR-L1 | RASQDVNTAVA | 159–169 | 163 165 166 167 | 161 162 164 168 169 |
| CDR-L2 | SASFLYS | 185–191 | 185 188 | 186 187 189 190 |
| CDR-L3 | QQHYTTPPT | 224–232 | 226 227 228 229 | 224 225 230 231 |

**Cross-validation: 22 of 24 zone-A residues sit inside a CDR.** Zones come
from 3D geometry, CDRs from sequence — the agreement is independent evidence
that both are right. The two outliers, **Y184 and R201**, are framework
residues that contact HER2 anyway (R201 buries 21 A^2). Worth a line in the
pitch: the paratope is not exactly the CDRs.

### Kabat numbering: recommend dropping the literature positions

The plan cites CDR-H3 affinity literature at "W98/F100c/A101/L102" (Kabat).
Our CDR-H3 is `WGGDGFYAMDY` = scFv 99–109. Kabat numbers this loop
95,96,97,98,99,100,100a,100b,100c,101,102 — so Kabat "98" is **D**, not W,
and the cited labels do not map cleanly onto this sequence.

**Do not translate those positions by hand.** Either run ANARCI
(`biomodals/modal_anarci.py`) to get an authoritative numbering, or — better
for today — ignore the literature positions entirely and take zone-A
proposals from the **AbMPNN x AbLang2 consensus**, which is generated in our
own numbering and needs no translation. The pitch is stronger for it:
structure-derived, not copied.

---

## 2026-10-31 — CORRECTION: zone A/B boundary was wrong

The hardcoded `ZONE_A` in the scan scripts held **24** residues — the
contact-**or**-burial list from `interface_residues.py`. The plan's definition
of zone A is **direct contact < 4.5 A**, which is **20** residues.

Four positions were therefore mislabelled A when they are interface rim (B):
**G101, D102, Y184, R201** — all bury surface on binding but make no contact
at 4.5 A.

Fixed in `annotate_positions.py`; `results/positions.csv` is now the single
source of truth, and `build_variants.py` reads zones from it rather than from
the stale `zone` columns in the scan CSVs. Consequences:

* **D102 is zone B**, consistent with the earlier liability analysis. The
  DG repair is a *rim* edit, not a paratope edit — lower risk than stated.
* Cross-validation improves: **20 of 20** zone-A residues are inside CDRs.
  Y184 and R201 are framework positions at the rim, burying 11 and 21 A^2
  without direct contact.

---

## 2026-10-31 — Phase 1c: AbMPNN (DONE)

`scripts/run_abmpnn.sh` -> 100 sequences in 46 s on CPU.
`scripts/parse_abmpnn.py` -> `results/abmpnn_scan.csv`
Score = `log2((n_mut + 0.5) / (n_wt + 0.5))` over 100 samples at T=0.2.
Native reference in the FASTA verified against the parent VH and VL
(output chain order is **A=VL then B=VH** — easy to transpose, don't).

Of 67 designed positions, AbMPNN varies only **23**. But where it varies, it
varies totally:

| position | parent kept | top proposal |
|---|---|---|
| Y227 | **0%** | Y227N (76/100) |
| N165 | **0%** | N165S (100/100) |
| F188 | **0%** | F188Y (94/100) |
| Y33 | **0%** | Y33W (100/100) |
| D102 | **0%** | D102Y (99/100) |
| T166 | **0%** | T166S (100/100) |
| T229 | **0%** | T229V (84/100) |
| W99 | 4% | W99E (84/100) |
| Y105 | 33% | Y105W (67/100) |

ProteinMPNN scores are negative log-likelihood, lower is better: the **native
scores 0.753, the designs score ~0.324**. AbMPNN thinks the real trastuzumab
paratope is a poor sequence for this backbone.

### The central finding of the day

**Two models, trained on completely different data, independently want to
dismantle the same paratope.** AbLang2 (natural repertoires, sequence) says
the CDRs are atypical. AbMPNN (antibody structures, backbone-conditioned)
says the residues don't suit the fold. Both are correct and both are useless
as instructions — trastuzumab's paratope is atypical **because** it was
affinity-matured and engineered against this epitope. Typicality is what was
optimised away.

Practical consequence: **neither model's zone-A proposals can be trusted on
its own score.** Only take zone-A edits where the two agree, and keep them to
the high-variance tier. This is also the pitch: we can show *why* naive
model-following fails on a matured antibody, with numbers.

### Consensus between AbLang2 and AbMPNN
Exact (position, mutation) agreement is rare — the useful signal:
**F188N** (zone A) and **D102S** (zone B, also the DG liability repair).
D102S is endorsed by both models *and* removes an isomerisation hotspot.

---

## 2026-10-31 — Phase 1b: ESM2 (DONE)

`scripts/esm2_scan.py` -> `results/esm2_scan.csv`. esm2_t33_650M, masked
marginals, 227 masked forward passes on CPU. **~10 min total: 8 min to
download the 2.6 GB checkpoint, ~2 min of inference.** Cached now, so a
re-run is ~2 min.

**302 positive substitutions** (A 105, B 119, C 78) — far more permissive
than AbLang2's 47, as expected from a general protein model.

### ESM2 is unreliable on antibody framework — demonstrated, not assumed

Top zone-C hit was **F218** (= VL F83), with **17 of 19** substitutions
scoring positive, best `F218E at +9.91`. AbLang2 scores the same mutation
**−9.15**. The structure settles it: **F218 is solvent-exposed, 84 A^2.**

A general PLM finds an exposed phenylalanine anomalous and wants it polar.
An antibody-trained model knows an exposed aromatic at that framework
position is canonical. **AbLang2 is right, ESM2 is wrong**, and this is
exactly the failure the multi-model design is meant to catch.

By contrast **A72R** is genuine consensus: ESM2 +5.18, AbLang2 +2.22,
exposed 44 A^2.

### Consequence: the framework has almost no defensible headroom

After the two-model veto, only **6** zone-C positions survive:
**28, 72, 79, 239, 240, 241**. And no 5 of them are mutually >6 A apart
(largest such set is 4: 28/72/239/241).

This is a real result, not a tooling failure: trastuzumab is already fully
humanised, so there is little germline distance left to close. Worth saying
in the pitch — most teams will assume framework mutations are free and
easy to find.

---

## 2026-10-31 — Phase 2: variant pool (DONE)

`scripts/build_variants.py` -> `results/candidate_mutations.csv`,
`results/variant_pool.csv`, `results/variant_pool.fasta`

### Candidate rule (two tiers of evidence)
* **zones A and B**: >= 2 of 3 sources must score it positive. Both models
  mis-call the matured paratope, so single-source zone-A calls are refused.
* **zone C**: one positive source suffices, provided **no source scores it
  below −2.0**. The veto is what excludes F218E. Without the relaxation the
  conservative tier is empty; without the veto it fills with ESM2 artefacts.

**58 candidate single mutations**: zone A 17, B 31, C 10, liability 6.

### Combination rules enforced
* >= 5 mutations, <= 260 aa (substitutions only, so always 242)
* at most one mutation per CDR loop
* at most one zone-A mutation per variant
* CB-CB > 10 A between chosen positions — relaxed to **6 A when both are
  zone C**, since distal surface framework positions share no environment
  (and 10 A makes the tier impossible)
* linker and the four conserved cysteines never touched

### Conservative tier needed redefinition
Pure zone C cannot reach 5 mutations (see above). The tier's *intent* is
"paratope untouched", so the pool was extended to zone-B **framework**
positions — no HER2 contact, outside the CDR loops. Available:
K30G/K30S, S97A, R201G.

### Pool: 53 variants, 0 rule violations

| tier | n | mutations each |
|---|---|---|
| conservative | 10 | 5 |
| rim | 17 | 6 |
| paratope | 10 | 6 (exactly one zone A) |
| liability | 7 | 5 (>=1 liability repair) |
| wildcard | 6 | 6 |
| wildcard-staple | 3 | 5 (G44C + Q234C + 3 zone C) |

All 53 are 242 aa with 5-6 mutations. Ready for Phase 3 cofolding.

---

## 2026-10-31 — Comparison with `SaniyaSani-patch-1`

That branch holds a **finished submission**: `final15.fasta`,
`final15_selection.csv`, `VALIDATION.txt`. They have already done our Phase 3.

### Their validation independently confirmed
15 unique ids, 15 unique sequences, all <= 242 aa, exactly one `(GGGGS)3`
linker each, 4 cysteines each (no new ones), canonical residues only, every
design >= 5 edits. **No rule failures.** Their stated edit distances (6, 5,
8, 9 for the SPR anchors) verified by pairwise alignment — a naive
position-wise diff reports 125/127 for the two 241 aa sequences because the
deletion shifts the frame; that was our measurement error, not their file.

### Convergence (independent, and strong)
| position | them | us |
|---|---|---|
| A72R | round-2 confirmed, +0.024 median paired delta | ESM2 +5.18 x AbLang2 +2.22 |
| N55 (NG deamidation) | mutated in all 11 computational | flagged zone B, free to fix |
| D102 (DG isomerisation) | mutated in all 11 | D102S dual-model consensus |

Two unrelated pipelines hit the same two liability motifs and the same
framework mutation.

### Divergence — concentrated risk in their pool
All 11 computational designs share the **same five positions**
(31, 55, 102, 163, 226); only the substituted residue varies at 31/55/102,
and two add A72R. **13 of 15 mutate >= 2 zone-A residues** — every
computational one changes both D163 (20 A^2, contacts HER2 598) and H226
(41 A^2, contacts 571/572). Our rules cap this at one zone-A edit.

If H226Y or D163K is wrong, all 11 fail together. This is the correlated
failure mode the organisers' "pick a diverse 15" is meant to prevent.

The 4 SPR anchors measure 0.63, 0.64, **0.87, 0.94** nM against a 0.82 nM
parent — two are *worse* than parent, and their own README notes the
platform is not comparable.

### Their methodology is better than ours in three ways
* **paired delta vs parent on the same seed** — cancels seed noise, which is
  stronger than comparing against a separate baseline run
* **pose_ok_rate** — fraction of seeds keeping a parent-like pose
* **epitope_jaccard** — overlap of the predicted epitope with the parent's.
  We had no equivalent. It directly catches a variant that scores well while
  docking somewhere else entirely. **Adopt this.**

### Correction to our own earlier advice on G103
`SPR_03` mutates **G103** (our 63 A^2 anchor) plus F104 and measures
**0.87 nM** vs the 0.82 nM parent — essentially unchanged. `SPR_04` mutates
W99 and F104 and measures 0.94 nM.

So "never touch G103" is too strong. Within a **coordinated full CDR-H3
redesign** the loop can be rebuilt around it and still bind sub-nM. The
warning holds for a *single point mutation* at G103 in an otherwise native
loop, which is a different proposition. State it that way in the pitch.

### Proposed merge
Run their validated scoring pipeline over our 53-variant pool, then build
the 15 from both sources: their top ~5 computational leaders plus our
conservative / liability / staple tiers, which cover positions they never
explored. Fixes their concentration and our lack of scores in one step.

---

## 2026-10-31 — Phase 3: Boltz screen

### Merged pool
`scripts/build_merged60.py` -> `results/merged60.{csv,fasta}`
60 variants: **22 ours, 11 Saniya's computational, 27 library**, chosen by
greedy max-min diversity over mutation sets so the three pools do not
duplicate each other. **44 distinct positions** covered; 27 variants touch no
zone-A residue.

Plus 3 controls: `parent`, `decoy_Y105A` (the 79 A^2 anchor deleted) and
`decoy_scrambleH3` (CDR-H3 shuffled). 63 sequences total.

Complex = 242 aa scFv (chain A) + 108 aa HER2 domain IV (chain B) = **350 aa**.

### Pilot (3 runs, seed 1) — the trim is validated

| | ipTM | pTM | pLDDT | d vs parent | epitope Jaccard |
|---|---|---|---|---|---|
| parent | 0.8907 | 0.9358 | 0.9125 | — | 1.00 |
| our_conservative_01 | 0.9150 | 0.9487 | 0.9213 | **+0.0243** | 0.895 |
| decoy_Y105A | 0.8502 | 0.9223 | 0.9224 | **−0.0405** | 0.75 |

**Parent cofolds at ipTM 0.89 against trimmed domain IV**, so the 490–607 trim
is sound and the variants can be interpreted against it.

**Noise band: +/- 0.04.** Deleting the single biggest anchor (Y105A, 79 A^2)
moves ipTM by only 0.0405. Everything between parent and "anchor deleted"
lives in that 0.04 window, so a variant within +/-0.04 of parent on one seed
is **not** distinguishable from it. Our +0.0243 variant is inside the band.

Y105A keeps **75% epitope overlap** — removing the anchor degrades the
interface without relocating it. So epitope Jaccard measures binding quality,
not mislocalisation, and a genuinely low Jaccard means something odd.

### Cost and timing (measured, not estimated)
~145 s per complex in steady state (first run 649 s including image build).
About $0.08 per run => **~$5 per seed, ~$10 for two**. Budget is $30; the
constraint is wall-clock, not money. Roughly 85% of each run is the GPU
container idling while it waits on the ColabFold MSA server.

### Gotcha: Modal app-create rate limit
Launching 6 `modal run` processes in the same second returns
`RESOURCE_EXHAUSTED: App create rate limit exceeded`. Failures are instant and
cost nothing, but 8 runs were lost before it was spotted. Fixed in
`run_boltz_screen.sh` with a 0-25 s jittered start and 4 retries at 30/60/90/
120 s **only** for rate-limit errors, so a genuine Boltz failure still fails
fast rather than being retried at $0.08 a go. Default concurrency reduced 6->4.

### Permission note
`ipsae.py` is cloned third-party code, so Claude cannot execute it under auto
mode — same as `biomodals` and AbMPNN. Run `bash scripts/run_ipsae.sh 1`
manually; it is local and free. `scripts/score_screen.py` is ours and runs
normally.

---

## Open questions / risks

1. **Kabat vs linear numbering.** CDR-H3 affinity literature cites W98/F100c/
   A101/L102 in Kabat. Our W99/F104/Y105 etc. are linear. Map explicitly
   before trusting any literature position, or we will mutate the wrong residue.
2. **AbLang2 and AbMPNN have no installed skill.** AbMPNN exists only inside
   the `germinal` skill's pipeline. Both need pulling from HuggingFace.
3. **GPU budget.** `boltz` skill quotes 100 complexes = 30–45 min = ~$8 on
   L40S. Phase 3 at 100 combos x 3 seeds = 300 complexes ≈ $25 of a $30 free
   tier, before Phase 0 baselines and the final full-ectodomain re-runs.
4. **Complex size.** 242 aa scFv + 630 aa full ECD = 872 aa. Trimmed to domain
   IV it is ~370 aa. Untrimmed runs risk CUDA OOM on L40S as well as cost.

---

## Results log

*(ipSAE / Boltz runs go here — nothing run yet)*

| run | construct | target | seeds | ipSAE_min mean | sd | notes |
|---|---|---|---|---|---|---|
| — | — | — | — | — | — | no GPU jobs run yet |
