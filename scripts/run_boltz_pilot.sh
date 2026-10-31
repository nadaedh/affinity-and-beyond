#!/usr/bin/env bash
# PILOT: 3 complexes, 1 seed. Measures real cost/time before committing the budget.
#
#   bash scripts/run_boltz_pilot.sh
#
# Afterwards check spend at https://modal.com/settings/usage and report the
# number -- that sets whether the full screen is 1 seed or 2.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT/biomodals"

SEED=1
OUT="$ROOT/results/boltz/seed${SEED}"
mkdir -p "$OUT"

# parent first (it is the control), then one real variant, then one decoy.
# The variant is whichever of ours the diversity picker actually selected.
VARIANT="$(basename "$(ls "$ROOT"/results/boltz_inputs/our_*.yaml | head -1)" .yaml)"
echo "pilot variant: $VARIANT"

for NAME in parent "$VARIANT" decoy_Y105A; do
    YAML="$ROOT/results/boltz_inputs/${NAME}.yaml"
    if [ ! -f "$YAML" ]; then
        echo "skip $NAME (no yaml)"; continue
    fi
    echo "=== $NAME (seed $SEED) ==="
    /usr/bin/time -p "$ROOT/.venv/bin/modal" run modal_boltz.py \
        --input-yaml "$YAML" \
        --run-name "${NAME}_s${SEED}" \
        --out-dir "$OUT" \
        --params-str "--use_msa_server --seed ${SEED}"
done

echo
echo "pilot done. outputs:"
find "$OUT" -name "*.cif" -o -name "*pae*" | head -20
