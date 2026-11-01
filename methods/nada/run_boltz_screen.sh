#!/usr/bin/env bash
# Full screen: every YAML in results/boltz_inputs at one seed, run in parallel.
#
#   bash scripts/run_boltz_screen.sh [SEED] [PARALLEL]
#   bash scripts/run_boltz_screen.sh 1 6
#
# Skips any run that already has a pae .npz, so it is safe to re-run after an
# interruption and safe to run straight after the pilot.
#
# Concurrency note: most of each run is spent waiting on the ColabFold MSA
# server, which rate-limits. 6 is a reasonable compromise; raising it may get
# you throttled rather than going faster.
set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SEED="${1:-1}"
PAR="${2:-4}"

OUT="$ROOT/results/boltz/seed${SEED}"
LOGS="$ROOT/results/boltz/logs/seed${SEED}"
mkdir -p "$OUT" "$LOGS"

run_one() {
    NAME="$1"; SEED="$2"; ROOT="$3"; OUT="$4"; LOGS="$5"
    DONE=$(find "$OUT/${NAME}_s${SEED}" -name "pae_*.npz" 2>/dev/null | head -1)
    if [ -n "$DONE" ]; then
        echo "skip  $NAME (already done)"
        return 0
    fi
    cd "$ROOT/biomodals"

    # Modal rate-limits app creation, so a simultaneous burst gets rejected.
    # Jitter the start, then retry with backoff on that specific failure.
    sleep $(( RANDOM % 25 ))

    for ATTEMPT in 1 2 3 4; do
        if "$ROOT/.venv/bin/modal" run modal_boltz.py \
                --input-yaml "$ROOT/results/boltz_inputs/${NAME}.yaml" \
                --run-name "${NAME}_s${SEED}" \
                --out-dir "$OUT" \
                --params-str "--use_msa_server --seed ${SEED}" \
                > "$LOGS/${NAME}.log" 2>&1; then
            echo "ok    $NAME (attempt $ATTEMPT)"
            return 0
        fi
        if grep -q "rate limit" "$LOGS/${NAME}.log"; then
            BACKOFF=$(( ATTEMPT * 30 + RANDOM % 20 ))
            echo "wait  $NAME rate-limited, retry in ${BACKOFF}s (attempt $ATTEMPT)"
            sleep "$BACKOFF"
        else
            echo "FAIL  $NAME  (see $LOGS/${NAME}.log)"
            return 1
        fi
    done
    echo "FAIL  $NAME  (rate-limited after 4 attempts)"
    return 1
}
export -f run_one

NAMES=$(cd "$ROOT/results/boltz_inputs" && ls *.yaml | sed 's/\.yaml$//')
TOTAL=$(echo "$NAMES" | wc -l | tr -d ' ')
echo "screening $TOTAL sequences at seed $SEED, $PAR at a time"
echo "start $(date +%H:%M:%S)"
echo

echo "$NAMES" | xargs -P "$PAR" -I {} bash -c 'run_one "$@"' _ {} "$SEED" "$ROOT" "$OUT" "$LOGS"

echo
echo "end   $(date +%H:%M:%S)"
NPZ=$(find "$OUT" -name "pae_*.npz" | wc -l | tr -d ' ')
echo "completed: $NPZ / $TOTAL"
FAILED=$(grep -L "pae" "$LOGS"/*.log 2>/dev/null | wc -l | tr -d ' ')
echo "check failures with: grep -il error $LOGS/*.log"
