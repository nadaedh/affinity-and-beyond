#!/usr/bin/env bash
# AbMPNN backbone-fixed redesign of the trastuzumab Fv against HER2 domain IV.
#
# Designs ONLY zone A + zone B positions (conserved cysteines excluded).
# HER2 (chain C) is held fixed, so proposals are epitope-aware.
#
# Run from the project root:  bash scripts/run_abmpnn.sh
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

PY="$ROOT/.venv/bin/python"
MPNN="$ROOT/tools/ProteinMPNN"
H="$MPNN/helper_scripts"
OUT="$ROOT/results/abmpnn"
TMP="$OUT/tmp"
mkdir -p "$TMP"

# design positions, written by scripts/prep_complex.py
VH_POS="$(sed -n 's/^B: //p' results/abmpnn_design_positions.txt)"
VL_POS="$(sed -n 's/^A: //p' results/abmpnn_design_positions.txt)"

# the complex must sit alone in a directory for parse_multiple_chains
mkdir -p "$TMP/pdbs"
cp data/structures/fv_her2_domainIV.pdb "$TMP/pdbs/"

"$PY" "$H/parse_multiple_chains.py" \
    --input_path "$TMP/pdbs" \
    --output_path "$TMP/parsed.jsonl"

# design chains B (VH) and A (VL); chain C (HER2) stays fixed
"$PY" "$H/assign_fixed_chains.py" \
    --input_path "$TMP/parsed.jsonl" \
    --output_path "$TMP/chains.jsonl" \
    --chain_list "B A"

# --specify_non_fixed: the positions listed are the ones that GET designed
"$PY" "$H/make_fixed_positions_dict.py" \
    --input_path "$TMP/parsed.jsonl" \
    --output_path "$TMP/positions.jsonl" \
    --chain_list "B A" \
    --position_list "$VH_POS, $VL_POS" \
    --specify_non_fixed

"$PY" "$MPNN/protein_mpnn_run.py" \
    --jsonl_path "$TMP/parsed.jsonl" \
    --chain_id_jsonl "$TMP/chains.jsonl" \
    --fixed_positions_jsonl "$TMP/positions.jsonl" \
    --path_to_model_weights "$MPNN/abmpnn_weights" \
    --model_name abmpnn \
    --out_folder "$OUT" \
    --num_seq_per_target 100 \
    --sampling_temp "0.2" \
    --seed 1 \
    --batch_size 1

echo
echo "done -> $OUT/seqs/"
ls -la "$OUT/seqs/" 2>/dev/null || true
