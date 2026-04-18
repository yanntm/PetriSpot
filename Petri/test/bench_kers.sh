#!/usr/bin/env bash
# bench_kers.sh — compare PNML vs KERS I/O performance on a large model.
#
# Runs petri64 on FamilyReunion in three modes:
#   1. Direct PNML -> P-flows  (baseline)
#   2. Direct PNML -> T-flows  (baseline)
#   3. Export PNML to KERS, then KERS -> P-flows
#   4. KERS -> T-flows
#
# loopLimit 500 is now the compiled-in default; no flag needed.
#
# Usage: bash bench_kers.sh [path/to/petri64] [path/to/model.pnml]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PETRI="${1:-$SCRIPT_DIR/../src/petri64}"
MODEL="${2:-$SCRIPT_DIR/../examples/FamilyReunion-PT-L00400M0040C020P020G001.pnml}"
KERS_FILE="/tmp/familyreunion.kers"
LOG="$SCRIPT_DIR/bench_kers_$(date +%Y%m%d_%H%M%S).log"

echo "=== KERS benchmark $(date) ===" | tee "$LOG"
echo "petri64 : $PETRI"               | tee -a "$LOG"
echo "model   : $MODEL"               | tee -a "$LOG"
echo                                  | tee -a "$LOG"

# --- File sizes ---
echo "--- File sizes ---"                                              | tee -a "$LOG"
ls -lh "$MODEL"                                                        | tee -a "$LOG"

echo ""                                                                | tee -a "$LOG"

# --- Step 1: Export PNML -> KERS ---
echo "--- Export PNML to KERS ---"                                    | tee -a "$LOG"
{ time "$PETRI" -i "$MODEL" --exportAsKERS="$KERS_FILE" -q ; } 2>&1  | tee -a "$LOG"
ls -lh "$KERS_FILE"                                                    | tee -a "$LOG"
echo ""                                                                | tee -a "$LOG"

# --- Step 2: PNML -> P-flows (baseline) ---
echo "--- PNML path: P-flows ---"                                     | tee -a "$LOG"
{ time "$PETRI" -i "$MODEL" --Pflows -q ; } 2>&1                     | tee -a "$LOG"
echo ""                                                                | tee -a "$LOG"

# --- Step 3: PNML -> T-flows (baseline) ---
echo "--- PNML path: T-flows ---"                                     | tee -a "$LOG"
{ time "$PETRI" -i "$MODEL" --Tflows -q ; } 2>&1                     | tee -a "$LOG"
echo ""                                                                | tee -a "$LOG"

# --- Step 4: KERS -> P-flows ---
echo "--- KERS path: P-flows ---"                                     | tee -a "$LOG"
{ time "$PETRI" --loadKERS="$KERS_FILE" --Pflows -q ; } 2>&1        | tee -a "$LOG"
echo ""                                                                | tee -a "$LOG"

# --- Step 5: KERS -> T-flows ---
echo "--- KERS path: T-flows ---"                                     | tee -a "$LOG"
{ time "$PETRI" --loadKERS="$KERS_FILE" --Tflows -q ; } 2>&1        | tee -a "$LOG"
echo ""                                                                | tee -a "$LOG"

echo "=== Done. Log written to $LOG ==="
