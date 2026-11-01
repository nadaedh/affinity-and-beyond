# Boltz prediction runs

Consolidated from 16 per-variant Modal console logs. The raw logs
contained no scores -- only progress output -- but they do record
something important:

**13 of 16 runs did not complete.** They failed with
`App create rate limit exceeded` from Modal, which rejects bursts of
simultaneous app creations. The scrambled-CDR-H3 negative control is
among the failures, so this pool has no decoy calibration.

The same failure hit nada's screen; the fix (jittered starts plus
backoff retries on that error only) is in
`methods/nada/run_boltz_screen.sh`.

| run | seed | completed | failure |
|---|---|---|---|
| CTRL_scrambled_H3 | - | **no** | Modal rate limit |
| V01_safe_base | 42 | yes |  |
| V02_base_D102S | - | **no** | Modal rate limit |
| V03_safe_alt | - | **no** | Modal rate limit |
| V04_ss_V168L | - | **no** | Modal rate limit |
| V05_ss_A106Y | - | **no** | Modal rate limit |
| V06_ss_combined | - | **no** | Modal rate limit |
| V07_deamid_ss | - | **no** | Modal rate limit |
| V08_ox_M107F | - | **no** | Modal rate limit |
| V09_ox_M107L | - | **no** | Modal rate limit |
| V10_interface_stack | - | **no** | Modal rate limit |
| V11_charge_D102R | - | **no** | Modal rate limit |
| V12_charge_double | - | **no** | Modal rate limit |
| V13_germline_test | - | **no** | Modal rate limit |
| V14_paratope_bet | 42 | yes |  |
| V15_combined_best | 42 | yes |  |

**3 completed: V01_safe_base, V14_paratope_bet, V15_combined_best**

Consequence: these 15 designs were selected largely from the
sequence-level scan in `scan_results.csv`, not from completed
structure predictions.
