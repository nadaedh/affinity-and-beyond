# affinity-and-beyond
Breast Cancer Awareness Hackathon: Protein design track

Trastuzumab scFv optimisation against HER2. Parent: 242 aa, VH–(GGGGS)3–VL,
0.82 nM KD. Variants need >=5 mutations from the parent, one chain, <=260 aa.

## Setup

Per-machine, not tracked in this repo.

**1. Protein design skills** (24 Claude Code skills from Adaptyv):

```
/plugin marketplace add adaptyvbio/protein-design-skills
/plugin install adaptyv@protein-design-skills
```

**2. Modal CLI** — serverless GPU, free tier covers a hackathon:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install modal
modal setup          # opens a browser, click Authorize
```

**3. biomodals** — the Modal run scripts the skills call:

```bash
git clone https://github.com/hgbrian/biomodals
```

Verify with `cd biomodals && modal run modal_boltzgen.py --help`.

Note: `modal token show` appears in the setup skill but no longer exists.
Use `modal profile current` to check auth.
