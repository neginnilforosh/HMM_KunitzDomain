#!/usr/bin/env python3
"""
STEP 4 — Build .class files from the current hmmsearch outputs
==============================================================

BUG FIXED:
    The original seq_set_N.class files were built from an old hmmsearch run
    and were never rebuilt after the sequence HMM was retrained.  As a result,
    the E-values in seq_set_1.class / seq_set_2.class did not match the
    current seq_pos_N.out / seq_neg_N.out files (0/106 positives matched).

    This script must be run AFTER step 3 to guarantee the .class files are
    always consistent with the HMM outputs that produced them.

OUTPUT FILES (written to Data/Processed/):
    seq_set_1.class      struct_set_1.class
    seq_set_2.class      struct_set_2.class

FORMAT (space-separated, one protein per line):
    <protein_id>  <label>  <evalue>
    label = 1 for positives, 0 for negatives
    evalue = full-sequence E-value from tblout; 10.0 if no hit (below threshold)

USAGE:
    python3 step4_build_class_files.py
    python3 step4_build_class_files.py --processed-dir /path/to/Data/Processed
"""

import argparse
import sys
from pathlib import Path


# ── Argument parsing ──────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--processed-dir", default="Data/Processed",
                   help="Path to the Data/Processed directory (default: Data/Processed)")
    return p.parse_args()


# ── File parsers ──────────────────────────────────────────────────────────────

def parse_tblout(filepath: Path) -> dict[str, float]:
    """
    Parse an hmmsearch --tblout file.
    Returns {protein_id: full_sequence_evalue} for all hits.
    Only the best-scoring line per protein is kept (tblout lists one line
    per target sequence with the best domain already summarised).
    """
    hits: dict[str, float] = {}
    with open(filepath) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.split()
            if len(cols) < 5:
                continue
            protein_id = cols[0]
            evalue = float(cols[4])          # full-sequence E-value (column 5)
            # Keep best (lowest) E-value if the same ID appears more than once
            if protein_id not in hits or evalue < hits[protein_id]:
                hits[protein_id] = evalue
    return hits


def read_fasta_ids(filepath: Path) -> list[str]:
    """Return ordered list of sequence IDs from a FASTA file."""
    ids = []
    with open(filepath) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    return ids


# ── Core builder ──────────────────────────────────────────────────────────────

def build_class_file(
    pos_fasta: Path,
    neg_fasta: Path,
    pos_out: Path,
    neg_out: Path,
    outfile: Path,
) -> None:
    """
    Build one .class file.

    Iterates over positives then negatives in FASTA order (preserving the
    same protein ordering used by the notebook's parse_class() function).
    Proteins with no hmmsearch hit receive E-value = 10.0 (above any
    biologically meaningful threshold, so they are classified as negatives
    at all thresholds we test).
    """
    pos_hits = parse_tblout(pos_out)
    neg_hits = parse_tblout(neg_out)

    pos_ids = read_fasta_ids(pos_fasta)
    neg_ids = read_fasta_ids(neg_fasta)

    lines = []

    for pid in pos_ids:
        ev = pos_hits.get(pid, 10.0)
        lines.append(f"{pid} 1 {ev}")

    for pid in neg_ids:
        ev = neg_hits.get(pid, 10.0)
        lines.append(f"{pid} 0 {ev}")

    with open(outfile, "w") as fh:
        fh.write("\n".join(lines) + "\n")

    n_pos = len(pos_ids)
    n_neg = len(neg_ids)
    n_pos_hits = sum(1 for pid in pos_ids if pid in pos_hits)
    n_neg_hits = sum(1 for pid in neg_ids if pid in neg_hits)
    print(f"  Written: {outfile}")
    print(f"    Positives: {n_pos} total, {n_pos_hits} with hits "
          f"({n_pos - n_pos_hits} assigned evalue=10.0)")
    print(f"    Negatives: {n_neg} total, {n_neg_hits} with hits "
          f"({n_neg - n_neg_hits} assigned evalue=10.0)")


# ── Validation ────────────────────────────────────────────────────────────────

def validate_class_file(class_file: Path, out_file: Path, label: str) -> bool:
    """
    Cross-check that every hit in the .out file appears with the same
    E-value in the .class file.  Reports mismatches.
    Returns True if all E-values match.
    """
    out_hits = parse_tblout(out_file)

    class_data: dict[str, float] = {}
    with open(class_file) as fh:
        for line in fh:
            parts = line.strip().split()
            if len(parts) >= 3:
                class_data[parts[0]] = float(parts[2])

    mismatches = 0
    for pid, ev in out_hits.items():
        if pid not in class_data:
            print(f"    MISSING in class file: {pid}")
            mismatches += 1
        elif class_data[pid] != ev:
            print(f"    MISMATCH {pid}: class={class_data[pid]:.2e}  out={ev:.2e}")
            mismatches += 1

    if mismatches == 0:
        print(f"    Validation OK — all {len(out_hits)} hits match {label}")
        return True
    else:
        print(f"    Validation FAILED — {mismatches} mismatches in {label}")
        return False


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()
    d = Path(args.processed_dir)

    if not d.is_dir():
        print(f"ERROR: processed directory not found: {d}")
        sys.exit(1)

    # Define the four (pos_fasta, neg_fasta, pos_out, neg_out, class_out) tuples
    builds = [
        (
            "seq_set_1",
            d / "pos_1.fasta",     d / "neg_1.fasta",
            d / "seq_pos_1.out",   d / "seq_neg_1.out",
            d / "seq_set_1.class",
        ),
        (
            "seq_set_2",
            d / "pos_2.fasta",     d / "neg_2.fasta",
            d / "seq_pos_2.out",   d / "seq_neg_2.out",
            d / "seq_set_2.class",
        ),
        (
            "struct_set_1",
            d / "pos_1.fasta",     d / "neg_1.fasta",
            d / "struct_pos_1.out", d / "struct_neg_1.out",
            d / "struct_set_1.class",
        ),
        (
            "struct_set_2",
            d / "pos_2.fasta",     d / "neg_2.fasta",
            d / "struct_pos_2.out", d / "struct_neg_2.out",
            d / "struct_set_2.class",
        ),
    ]

    print("=== STEP 4: Build .class files from current hmmsearch outputs ===")
    print()

    # Check all required input files exist before writing any output
    missing = []
    for name, pf, nf, po, no, _ in builds:
        for fpath in (pf, nf, po, no):
            if not fpath.exists():
                missing.append(str(fpath))
    if missing:
        print("ERROR: the following required files are missing:")
        for m in missing:
            print(f"  {m}")
        print("Run step 3 (step3_run_hmmsearch.sh) first.")
        sys.exit(1)

    # Build all four class files
    all_valid = True
    for name, pos_fasta, neg_fasta, pos_out, neg_out, outfile in builds:
        print(f"Building {name}.class ...")
        build_class_file(pos_fasta, neg_fasta, pos_out, neg_out, outfile)

        # Validate both pos and neg portions
        pos_label = outfile.stem + " (positives)"
        neg_label = outfile.stem + " (negatives)"
        ok_pos = validate_class_file(outfile, pos_out, pos_label)
        ok_neg = validate_class_file(outfile, neg_out, neg_label)
        if not (ok_pos and ok_neg):
            all_valid = False
        print()

    # Final summary
    print("=== Summary ===")
    for _, _, _, _, _, outfile in builds:
        n_lines = sum(1 for _ in open(outfile))
        n_pos = sum(1 for l in open(outfile) if l.split()[1] == "1")
        n_neg = sum(1 for l in open(outfile) if l.split()[1] == "0")
        print(f"  {outfile.name:<25}  {n_lines} entries  "
              f"(pos={n_pos}, neg={n_neg})")

    print()
    if all_valid:
        print("All .class files built and validated successfully.")
        print("Proceed to the MCC_vs_Evalue_Plot.ipynb notebook.")
    else:
        print("WARNING: validation failures detected — check output above.")
        sys.exit(1)


if __name__ == "__main__":
    main()
