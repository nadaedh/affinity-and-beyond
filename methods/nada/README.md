# nada

Structure-first: map the interface, point three models at the zones each is
competent in, make them veto each other, then screen with decoy controls.

**Start here:** [`lab-notebook.md`](lab-notebook.md) — the full decision log,
including the mistakes and corrections.

| stage | scripts |
|---|---|
| interface map | `interface_residues.py`, `zones_and_liabilities.py`, `annotate_positions.py` |
| sequence scans | `ablang2_scan.py`, `esm2_scan.py`, `run_abmpnn.sh` + `parse_abmpnn.py` |
| variant assembly | `build_variants.py`, `prep_complex.py` |
| screening | `make_boltz_yamls.py`, `run_boltz_pilot.sh`, `run_boltz_screen.sh` |
| scoring | `run_ipsae.sh`, `score_screen.py` |
| selection | `select_final15.py`, `build_merged60.py`, `compare_pools.py` |

| output | what |
|---|---|
| `merged60.csv` | the 60 screened variants, pooled from three team sources |
| `nada_pool_final15.fasta` | this pool's 15, spanning 32 positions — **not submitted** |

## Two things this pipeline contributed

**Decoy calibration.** The only pool that scored deliberately broken controls
(Y105A and scrambled CDR-H3) alongside the designs, which is what made the
noise floor measurable.

**A rate-limit fix.** `run_boltz_screen.sh` jitters starts and retries with
backoff *only* on Modal's app-create rate limit, so a genuine failure still
fails fast. The same error silently killed 13 of 16 runs in another pool.

## Limitation

Single seed per variant. Individually underpowered — see
[`../../analysis/README.md`](../../analysis/README.md) section 3.
