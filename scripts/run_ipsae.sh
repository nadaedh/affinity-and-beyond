#!/usr/bin/env bash
# Run ipSAE over every finished Boltz prediction.
#
#   bash scripts/run_ipsae.sh [SEED]
#
# ipsae.py writes its output files next to the input .cif. Cheap and local --
# no GPU, no Modal spend.
set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SEED="${1:-1}"
PY="$ROOT/.venv/bin/python"
IPSAE="$ROOT/tools/IPSAE/ipsae.py"

n=0
for PAE in $(find "$ROOT/results/boltz/seed${SEED}" -name "pae_*.npz" | sort); do
    DIR="$(dirname "$PAE")"
    CIF="${DIR}/$(basename "$PAE" | sed 's/^pae_//; s/\.npz$/.cif/')"
    [ -f "$CIF" ] || { echo "no cif for $PAE"; continue; }
    # ipsae.py writes beside the cif; skip if already scored
    if ls "${DIR}"/*ipsae*.txt >/dev/null 2>&1; then continue; fi
    ( cd "$DIR" && "$PY" "$IPSAE" "$(basename "$PAE")" "$(basename "$CIF")" 10 10 >/dev/null 2>&1 )
    n=$((n+1))
done
echo "ran ipsae on $n predictions"
echo "output files:"
find "$ROOT/results/boltz/seed${SEED}" -iname "*ipsae*" | head -5
