#!/bin/bash
# =============================================================================
# STEP 2 — Build the structure-based HMM from the PDB structural alignment
# =============================================================================
#
# REQUIRES:  hmmer >= 3.1  (hmmbuild)
#
# USAGE:  bash step2_build_struct_hmm.sh
# =============================================================================

set -euo pipefail

PROCESSED="Data/Processed"
STRUCT_ALN="$PROCESSED/pdb_kunitz_nr_clean.ali"
STRUCT_HMM="$PROCESSED/structure_based.hmm"

echo "=== STEP 2: Build structure-based HMM from PDB structural alignment ==="
echo ""

# ── 2a. Verify input alignment ────────────────────────────────────────────────
if [ ! -f "$STRUCT_ALN" ]; then
    echo "ERROR: $STRUCT_ALN not found."
    exit 1
fi

N_SEQS=$(grep -c "^>" "$STRUCT_ALN")
echo "Structural alignment: $STRUCT_ALN  ($N_SEQS sequences)"

if [ "$N_SEQS" -gt 35 ]; then
    echo "WARNING: alignment has $N_SEQS sequences — expected ~24."
fi

# Check IDs are PDB-format (should NOT start with >sp|)
FIRST_ID=$(grep "^>" "$STRUCT_ALN" | head -1)
if echo "$FIRST_ID" | grep -q "^>sp|"; then
    echo "ERROR: alignment contains UniProt sequences (IDs start with >sp|)."
    echo "       pdb_kunitz_nr_clean.ali should contain PDB chain IDs only."
    exit 1
fi

# ── 2b. Build HMM ─────────────────────────────────────────────────────────────
echo ""
echo "Building structure HMM with hmmbuild..."

# Use -n (short flag) — older macOS hmmbuild does not accept --name (long form)
hmmbuild -n structure_based "$STRUCT_HMM" "$STRUCT_ALN"

if [ ! -s "$STRUCT_HMM" ]; then
    echo "ERROR: hmmbuild produced no output."
    exit 1
fi

echo "HMM written: $STRUCT_HMM"

# ── 2c. Sanity check ──────────────────────────────────────────────────────────
NSEQ=$(grep "^NSEQ" "$STRUCT_HMM" | awk '{print $2}')
LENG=$(grep "^LENG" "$STRUCT_HMM" | awk '{print $2}')
echo ""
echo "HMM summary:"
echo "  NSEQ (training sequences): $NSEQ  (expected ~24)"
echo "  LENG (model length):       $LENG  (expected ~55–65 aa)"

echo ""
echo "Step 2 complete. Proceed to: bash step3_run_hmmsearch.sh"
