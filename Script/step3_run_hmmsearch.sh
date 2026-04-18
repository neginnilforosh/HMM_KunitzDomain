#!/bin/bash
# =============================================================================
# STEP 3 — Run hmmsearch for both HMMs against all test FASTA sets
# =============================================================================
#
# Produces 8 .out files (tblout format):
#
#   seq_pos_1.out   struct_pos_1.out   — hmmsearch against pos_1.fasta
#   seq_pos_2.out   struct_pos_2.out   — hmmsearch against pos_2.fasta
#   seq_neg_1.out   struct_neg_1.out   — hmmsearch against neg_1.fasta
#   seq_neg_2.out   struct_neg_2.out   — hmmsearch against neg_2.fasta
#
# IMPORTANT — why we search pos and neg SEPARATELY:
#   fix_class_files.py (step 4) needs to know which FASTA set a protein came
#   from so it can assign the correct label (1 for pos, 0 for neg).  If we
#   merged pos and neg into one FASTA and searched once, we could not cleanly
#   assign labels afterward.
#
# REQUIRES:  hmmer >= 3.3  (hmmsearch)
#            Both HMMs from steps 1 and 2 must be present.
#
# USAGE:  bash step3_run_hmmsearch.sh
#         Optional: bash step3_run_hmmsearch.sh --threads 8
# =============================================================================

set -euo pipefail

PROCESSED="Data/Processed"
SEQ_HMM="$PROCESSED/sequence_based.hmm"
STRUCT_HMM="$PROCESSED/structure_based.hmm"

# Test FASTA files (produced by the cross-validation split upstream)
POS1="$PROCESSED/pos_1.fasta"
POS2="$PROCESSED/pos_2.fasta"
NEG1="$PROCESSED/neg_1.fasta"
NEG2="$PROCESSED/neg_2.fasta"

# Number of CPU threads (override with --threads N)
THREADS=4
if [[ "${1:-}" == "--threads" && -n "${2:-}" ]]; then
    THREADS="$2"
fi

echo "=== STEP 3: Run hmmsearch for both HMMs against all test sets ==="
echo "  Threads: $THREADS"
echo ""

# ── Helper: check file exists and is non-empty ────────────────────────────────
require_file() {
    if [ ! -s "$1" ]; then
        echo "ERROR: required file not found or empty: $1"
        exit 1
    fi
}

require_file "$SEQ_HMM"
require_file "$STRUCT_HMM"
require_file "$POS1"
require_file "$POS2"
require_file "$NEG1"
require_file "$NEG2"

# ── Sanity: verify HMMs were built from seeds, not the full dataset ───────────
SEQ_NSEQ=$(grep "^NSEQ" "$SEQ_HMM" | awk '{print $2}')
if [ "$SEQ_NSEQ" -gt 30 ]; then
    echo "ERROR: sequence_based.hmm has NSEQ=$SEQ_NSEQ."
    echo "       It must be built from the 24 PDB seeds (step 1), not from the"
    echo "       211-sequence ok_kunitz_clean set.  Re-run step 1."
    exit 1
fi
echo "Sequence HMM: NSEQ=$SEQ_NSEQ (OK)"
STRUCT_NSEQ=$(grep "^NSEQ" "$STRUCT_HMM" | awk '{print $2}')
echo "Structure HMM: NSEQ=$STRUCT_NSEQ (OK)"
echo ""

# ── Helper: run one hmmsearch ─────────────────────────────────────────────────
run_search() {
    local hmm="$1"
    local fasta="$2"
    local outfile="$3"
    local label="$4"

    echo "  Running: $label"
    echo "    HMM:   $hmm"
    echo "    FASTA: $fasta  ($(grep -c "^>" "$fasta") sequences)"
    echo "    OUT:   $outfile"

    hmmsearch \
        --tblout "$outfile" \
        --noali \
        --cpu "$THREADS" \
        -E 10 \
        "$hmm" \
        "$fasta" \
        > /dev/null

    local nhits
    nhits=$(grep -v "^#" "$outfile" | wc -l | tr -d ' ')
    echo "    Hits (E<=10): $nhits"
    echo ""
}

# ── Sequence-based searches ───────────────────────────────────────────────────
echo "--- Sequence-based HMM ---"
run_search "$SEQ_HMM" "$POS1" "$PROCESSED/seq_pos_1.out" "seq vs pos_1"
run_search "$SEQ_HMM" "$POS2" "$PROCESSED/seq_pos_2.out" "seq vs pos_2"
run_search "$SEQ_HMM" "$NEG1" "$PROCESSED/seq_neg_1.out" "seq vs neg_1"
run_search "$SEQ_HMM" "$NEG2" "$PROCESSED/seq_neg_2.out" "seq vs neg_2"

# ── Structure-based searches ──────────────────────────────────────────────────
echo "--- Structure-based HMM ---"
run_search "$STRUCT_HMM" "$POS1" "$PROCESSED/struct_pos_1.out" "struct vs pos_1"
run_search "$STRUCT_HMM" "$POS2" "$PROCESSED/struct_pos_2.out" "struct vs pos_2"
run_search "$STRUCT_HMM" "$NEG1" "$PROCESSED/struct_neg_1.out" "struct vs neg_1"
run_search "$STRUCT_HMM" "$NEG2" "$PROCESSED/struct_neg_2.out" "struct vs neg_2"

# ── Summary ───────────────────────────────────────────────────────────────────
echo "=== Step 3 complete. Output files: ==="
for f in seq_pos_1 seq_pos_2 seq_neg_1 seq_neg_2 \
          struct_pos_1 struct_pos_2 struct_neg_1 struct_neg_2; do
    FPATH="$PROCESSED/${f}.out"
    NHITS=$(grep -v "^#" "$FPATH" | wc -l | tr -d ' ')
    echo "  $FPATH  ($NHITS hits)"
done
echo ""
echo "Proceed to step 4 (fix_class_files.py) to build .class files."
