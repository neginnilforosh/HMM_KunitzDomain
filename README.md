# 🧬 Profile HMM Detection of Kunitz-Type Protease Inhibitor Domains (PF00014)

[![Python](https://img.shields.io/badge/python-3.10+-blue)]()
[![HMMER](https://img.shields.io/badge/HMMER-3.4-green)]()
[![Platform](https://img.shields.io/badge/platform-linux%20%7C%20macOS-lightgrey)]()

> Building and comparing **structure-based** vs **sequence-based** Profile Hidden Markov Models for Kunitz/BPTI-type protease inhibitor domain (Pfam PF00014) detection, evaluated by 2-fold cross-validation on Swiss-Prot data.

---

## 📋 Table of Contents

- [Overview](#-overview)
- [Repository Structure](#-repository-structure)
- [Methodology](#-methodology)
- [Results](#-results)
- [Installation](#️-installation)
- [Reproducing the Pipeline](#-reproducing-the-pipeline)
- [Data Availability](#-data-availability)

---

## 📌 Overview

This project constructs and compares two Profile HMMs for detecting the Kunitz/BPTI-type protease inhibitor domain (Pfam PF00014):

| Model               | Alignment Source                          | Seeds                                          | Tool                  |
|---------------------|-------------------------------------------|------------------------------------------------|-----------------------|
| **Structure-based** | PDBeFold structural superposition (SSM)   | 24 PDB chains (Pfam PF00014, ≤3.5 Å, 45–80 aa) | PDBeFold → `hmmbuild` |
| **Sequence-based**  | MUSCLE multiple sequence alignment        | 206 Swiss-Prot Kunitz domain sequences         | MUSCLE → `hmmbuild`   |

Both models are evaluated on the same Swiss-Prot positive/negative datasets using **2-fold cross-validation**, with MCC as the primary metric.

### Key Finding

Both models achieve near-perfect classification (peak MCC ≈ 0.995), but the **structure-based HMM reaches zero false positives at a 1,000× less strict threshold** (1e-5 vs 1e-8 in Fold 1). The structure-based model also assigns substantially more extreme E-values to true Kunitz sequences (best 8.6e-194 vs 7.5e-101), confirming that incorporating 3D structural information produces more discriminative emission probabilities.

---

## 📁 Repository Structure

```
HMM_KunitzDomain/
│
├── Data/
│   ├── Raw/
│   │   ├── pdb_ids.txt                    # 130 PDB IDs used (Pfam PF00014, ≤3.5Å, 45–80 aa)
│   │   └── kunitz_pdb_nr90.fasta          # 35 chains after CD-HIT 90% clustering
│   │
│   └── Processed/
│       ├── pdb_kunitz_nr_clean.ali        # ★ Structure-based alignment (PDBeFold SSM, 24 chains)
│       ├── seq_domains_aln.ali            # ★ Sequence-based alignment (MUSCLE, 206 domains)
│       ├── structure_based.hmm            # ★ Structure-based Profile HMM (LENG=58, NSEQ=24, EFFN=4.09)
│       ├── sequence_based.hmm             # ★ Sequence-based Profile HMM (LENG=58, NSEQ=206, EFFN=14.74)
│       ├── kunitz_true_seeds.fasta        # 24 PDBeFold seed domain sequences
│       ├── kunitz_domains_swssprot.fasta  # 206 Swiss-Prot Kunitz domain subsequences
│       ├── all_kunitz.fasta               # Merged human + non-human Kunitz positives
│       ├── pos_1.fasta, pos_2.fasta       # Positive cross-validation folds (211 total)
│       ├── neg_1.fasta, neg_2.fasta       # Negative cross-validation folds (~287K total)
│       └── contaminated_ids.txt           # 8 Swiss-Prot IDs ≥95% identical to PDB seeds
│
├── Scripts/
│   ├── performance.py                     # Confusion matrix, MCC, TPR, PPV, FPR computation
│   ├── build_class_files.py               # Merge hmmsearch outputs → .class files
│   ├── extract_kunitz_domains.py          # Extract Kunitz domain coordinates from hmmsearch domtbl
│   ├── MCC_vs_Evalue_Plot.ipynb           # Figure: MCC curves for both models across thresholds
│   ├── Confusion_Matrix_Table.ipynb       # Figure: Confusion matrices at optimal threshold
│   └── Evalue_Distribution_Plot.ipynb     # Figure: E-value score distributions
│
├── Results/
│   ├── all_performance_results.txt        # Full MCC/TPR/PPV/FPR sweep (all thresholds, both folds)
│   └── Figures/
│       ├── Figure_MCC_vs_Evalue.png       # MCC curves: structure-based vs sequence-based
│       ├── Figure_Confusion_Matrices.png  # Confusion matrices at optimal threshold
│       ├── Figure_Evalue_Distributions.png    # Score separation histograms
│       └── Figure_Evalue_Folds_Comparison.png # Fold consistency comparison
│
├── environment.yml                        # Conda environment specification
├── .gitignore
└── README.md
```

Files marked ★ are the core scientific deliverables — the alignments and trained HMMs.

---

## 🔬 Methodology

### 1. Data Preparation

**Positive set** — UniProtKB/Swiss-Prot release 2025_03:

- Candidates extracted under three annotation sources: Pfam **PF00014**, InterPro **IPR036880**, and Prosite **PS00280** / **PS50279**
- Queried separately for human (`organism_id:9606`) and non-human, then merged into `all_kunitz.fasta`
- **211 sequences** retained after CD-HIT clustering at 90% identity
- A BLASTp filter identified **8 sequences** with ≥95% identity to the PDB structural seeds
- The 8 contaminated sequences were **excluded from the structural training alignment** to prevent data leakage, but **retained in the cross-validation folds**
- Final: **211 positives** split into 2 folds (`pos_1.fasta`, `pos_2.fasta`)

**Negative set** — All remaining Swiss-Prot entries lacking Pfam PF00014 annotation (≈287,000 sequences), split into `neg_1.fasta`, `neg_2.fasta`.

> ⚠️ **Note:** The negative set is filtered only on PF00014. Some sequences carry Kunitz annotations under InterPro IPR036880 or Prosite PS00280/PS50279 but not Pfam — these land in the negative set and are responsible for most of the residual "false positives" at permissive thresholds (see *Discussion* in the report).

**PDB structures** — RCSB Advanced Search with three filters applied simultaneously: Pfam PF00014, resolution ≤3.5 Å, and chain length 45–80 aa:

- 130 PDB entries downloaded (`pdb_ids.txt`)
- CD-HIT clustering at 90% identity → 35 representative chains (`kunitz_pdb_nr90.fasta`)
- PDBeFold structural superposition (SSM) + visual inspection of the canonical Kunitz fold → **24 chains** retained for the structural alignment

### 2. Structure-Based HMM

```bash
# Structural alignment performed via PDBeFold (SSM, Q-score 0.7, all-vs-all)
# → pdb_kunitz_nr_clean.ali (24 sequences)
hmmbuild Data/Processed/structure_based.hmm Data/Processed/pdb_kunitz_nr_clean.ali
# → LENG = 58, NSEQ = 24, EFFN = 4.09
```

### 3. Sequence-Based HMM

```bash
# Step 3a: locate domain coordinates using the structure-based HMM
hmmsearch --domtblout Data/Processed/domain_hits.tbl --max -Z 1000 \
          Data/Processed/structure_based.hmm Data/Processed/ok_kunitz_clean.fasta

# Step 3b: extract domain subsequences (40–100 aa) from full-length Swiss-Prot positives
python3 Scripts/extract_kunitz_domains.py
# → kunitz_domains_swssprot.fasta (206 domains)

# Step 3c: align with MUSCLE v3.8 (ClustalW format) and build the HMM
muscle -in Data/Processed/kunitz_domains_swssprot.fasta \
       -out Data/Processed/seq_domains_aln.ali -clwstrict
hmmbuild Data/Processed/sequence_based.hmm Data/Processed/seq_domains_aln.ali
# → LENG = 58, NSEQ = 206, EFFN = 14.74
```

> The structure-based HMM is used to identify domain boundaries so that only the Kunitz region (not full-length Swiss-Prot sequences) feeds the sequence-based alignment. This keeps both models operating on equivalent 58-position windows.

### 4. HMM Search (2-fold cross-validation)

```bash
# 8 searches total: 2 models × 2 folds × (positives + negatives)
for fold in 1 2; do
    hmmsearch -Z 1000 --max --tblout Data/Processed/struct_pos_${fold}.out \
              Data/Processed/structure_based.hmm Data/Processed/pos_${fold}.fasta
    hmmsearch -Z 1000 --max --tblout Data/Processed/struct_neg_${fold}.out \
              Data/Processed/structure_based.hmm Data/Processed/neg_${fold}.fasta
    hmmsearch -Z 1000 --max --tblout Data/Processed/seq_pos_${fold}.out \
              Data/Processed/sequence_based.hmm Data/Processed/pos_${fold}.fasta
    hmmsearch -Z 1000 --max --tblout Data/Processed/seq_neg_${fold}.out \
              Data/Processed/sequence_based.hmm Data/Processed/neg_${fold}.fasta
done
```

> **Note on `-Z 1000`:** This sets the effective database size used for E-value computation. Because the actual per-fold database is ≈287,000 sequences, the reported E-values are deflated by a factor of ≈287 relative to a correctly calibrated search. MCC rankings and confusion-matrix counts are unaffected (E-value ordering is preserved), but the absolute thresholds reported here would shift roughly 2–3 orders of magnitude toward more permissive values under proper calibration. They are therefore not directly comparable to Pfam's standard reporting thresholds.

### 5. Performance Evaluation

```bash
# Build .class files (format: ID label E-value).
# Sequences absent from hmmsearch output are assigned a placeholder E-value of 10.
python3 Scripts/build_class_files.py

# Sweep E-value thresholds 10⁻¹ → 10⁻¹²
for model in struct seq; do
    for fold in 1 2; do
        for i in $(seq 1 12); do
            python3 Scripts/performance.py \
                Data/Processed/${model}_set_${fold}.class 1e-$i
        done
    done
done
```

The Matthews Correlation Coefficient (MCC) is the primary metric, since accuracy alone is uninformative under the strong class imbalance (~287,000 negatives vs. ~105 positives per fold).

---

## 📊 Results

### MCC at Key Thresholds

| Threshold | Struct Fold 1 | Struct Fold 2 | Seq Fold 1   | Seq Fold 2  |
|-----------|---------------|---------------|--------------|-------------|
| 1e-3      | 0.9768        | 0.9812        | 0.9643       | 0.9726      |
| 1e-4      | 0.9859        | 0.9714        | 0.9816       | 0.9770      |
| **1e-5**  | **0.9953** ⭐  | 0.9808       | 0.9859       | 0.9906      |
| **1e-6**  | 0.9905        | **0.9856** ⭐ | 0.9906       | **0.9953** ⭐|
| 1e-7      | 0.9905        | 0.9808        | 0.9905       | 0.9952      |
| **1e-8**  | 0.9905        | 0.9808        | **0.9953** ⭐ | 0.9952      |
| 1e-10     | 0.9809        | 0.9710        | 0.9905       | 0.9904      |
| 1e-12     | 0.9761        | 0.9611        | 0.9857       | 0.9808      |

### Classification Metrics at Optimal Thresholds

|         | Struct F1 (1e-5) | Struct F2 (1e-6) | Seq F1 (1e-8) | Seq F2 (1e-6) |
|---------|------------------|------------------|---------------|---------------|
| **TN**  | 287,115          | 287,114          | 287,115       | 287,113       |
| **FP**  | 0                | 0                | 0             | 1             |
| **FN**  | 1                | 3                | 1             | 0             |
| **TP**  | 105              | 102              | 105           | 105           |
| **MCC** | 0.9953           | 0.9856           | 0.9953        | 0.9953        |
| **TPR** | 0.9906           | 0.9714           | 0.9906        | 1.0000        |
| **PPV** | 1.0000           | 1.0000           | 1.0000        | 0.9906        |

### Key Figures

| Figure | Description |
|---|---|
| [Figure_MCC_vs_Evalue.png](Results/Figures/Figure_MCC_vs_Evalue.png) | MCC curves for both models across E-value thresholds |
| [Figure_Confusion_Matrices.png](Results/Figures/Figure_Confusion_Matrices.png) | Confusion matrices at optimal threshold |
| [Figure_Evalue_Distributions.png](Results/Figures/Figure_Evalue_Distributions.png) | Score separation histograms (positives vs. negatives) |

### Scientific Conclusion

Both models achieve near-perfect classification over more than **574,000** evaluated sequences, but the **structure-based HMM demonstrates superior discriminative power**: it achieves zero false positives at threshold 1e-5, while the sequence-based model requires the **1,000× stricter** threshold of 1e-8 to do the same in Fold 1 (both models share the 1e-6 optimum in Fold 2). The structure-based model also assigns more extreme E-values to true Kunitz sequences (best 8.6e-194 vs. 7.5e-101 — nearly 93 orders of magnitude tighter).

The few residual errors trace mostly to **annotation gaps** in Swiss-Prot rather than model failures: false positives at permissive thresholds are typically proteins (e.g., APP, TFPI) carrying Kunitz annotations under InterPro or Prosite but not Pfam PF00014. Using three annotation sources to build the positive set reduces this problem but does not eliminate it.

---

## ⚙️ Installation

```bash
git clone https://github.com/neginnilforosh/HMM_KunitzDomain.git
cd HMM_KunitzDomain
conda env create -f environment.yml
conda activate hmm-kunitz
```

**External tools required:**

- HMMER **v3.4** (`hmmbuild`, `hmmsearch`)
- MUSCLE **v3.8** (ClustalW output mode)
- CD-HIT **v4.8**
- BLAST+ **v2.12** (`makeblastdb`, `blastp`)
- Python **3.10+** with BioPython, NumPy, Matplotlib

---

## 🔁 Reproducing the Pipeline

### Step 1 — Download raw data

```bash
# Positive set from UniProtKB/Swiss-Prot (release 2025_03)
# Human:     https://www.uniprot.org/uniprotkb?query=(database:pfam+PF00014)+AND+reviewed:true+AND+organism_id:9606
# Non-human: https://www.uniprot.org/uniprotkb?query=(database:pfam+PF00014)+AND+reviewed:true+NOT+organism_id:9606
# Negative:  https://www.uniprot.org/uniprotkb?query=reviewed:true+NOT+(database:pfam+PF00014)
# → Save as: Data/Raw/human_kunitz.fasta, Data/Raw/nothuman_kunitz.fasta, Data/Raw/swissprot_notkunitz.fasta

# PDB structures — use pdb_ids.txt with RCSB batch download
# https://www.rcsb.org/downloads → paste contents of Data/Raw/pdb_ids.txt → select PDB format
# → Save all .pdb files to: Data/Raw/pdb_structures/
```

### Step 2 — Process and cluster the positive set

```bash
# Combine human and non-human Kunitz sequences
cat Data/Raw/human_kunitz.fasta Data/Raw/nothuman_kunitz.fasta > Data/Processed/all_kunitz.fasta

# Cluster at 90% identity with CD-HIT → 211 sequences
cd-hit -i Data/Processed/all_kunitz.fasta -o Data/Processed/ok_kunitz.fasta -c 0.90 -n 5

# BLASTp contamination check — flag sequences ≥95% identical to PDB seeds
makeblastdb -in Data/Processed/kunitz_true_seeds.fasta -dbtype prot -out Data/Raw/seeds_db
blastp -query Data/Processed/ok_kunitz.fasta -db Data/Raw/seeds_db \
       -outfmt 6 -out Data/Raw/blast_results.txt -evalue 1e-3
awk '$3 >= 95 {print $1}' Data/Raw/blast_results.txt | sort -u > Data/Processed/contaminated_ids.txt
# → 8 sequences flagged

# IMPORTANT: the 8 contaminated sequences are excluded ONLY from the structural training
# alignment (to avoid data leakage). They are RETAINED in the cross-validation folds,
# so the full positive set used for evaluation contains all 211 sequences.

# Split positives into 2 folds (all 211 sequences kept)
python3 -c "
from Bio import SeqIO
records = list(SeqIO.parse('Data/Processed/ok_kunitz.fasta', 'fasta'))
SeqIO.write(records[::2],  'Data/Processed/pos_1.fasta', 'fasta')
SeqIO.write(records[1::2], 'Data/Processed/pos_2.fasta', 'fasta')
print(f'pos_1: {len(records[::2])} | pos_2: {len(records[1::2])}')
"

# Split negatives into 2 folds
python3 -c "
from Bio import SeqIO
records = list(SeqIO.parse('Data/Raw/swissprot_notkunitz.fasta', 'fasta'))
SeqIO.write(records[::2],  'Data/Processed/neg_1.fasta', 'fasta')
SeqIO.write(records[1::2], 'Data/Processed/neg_2.fasta', 'fasta')
print(f'neg_1: {len(records[::2])} | neg_2: {len(records[1::2])}')
"
```

### Step 3 — Build the HMMs

```bash
# Structure-based HMM — alignment already provided in the repo
hmmbuild Data/Processed/structure_based.hmm Data/Processed/pdb_kunitz_nr_clean.ali

# Sequence-based HMM — extract domain regions from Swiss-Prot positives
# Step 3a: locate domain coordinates using the structure-based HMM
hmmsearch --domtblout Data/Processed/domain_hits.tbl --max -Z 1000 \
          Data/Processed/structure_based.hmm Data/Processed/ok_kunitz.fasta

# Step 3b: extract 40–100 aa domain subsequences (reads domain_hits.tbl automatically)
python3 Scripts/extract_kunitz_domains.py

# Step 3c: align with MUSCLE and build the HMM
muscle -in Data/Processed/kunitz_domains_swssprot.fasta \
       -out Data/Processed/seq_domains_aln.ali -clwstrict
hmmbuild Data/Processed/sequence_based.hmm Data/Processed/seq_domains_aln.ali
```

### Step 4 — Run hmmsearch (8 searches total)

```bash
for fold in 1 2; do
    hmmsearch -Z 1000 --max --tblout Data/Processed/struct_pos_${fold}.out \
              Data/Processed/structure_based.hmm Data/Processed/pos_${fold}.fasta
    hmmsearch -Z 1000 --max --tblout Data/Processed/struct_neg_${fold}.out \
              Data/Processed/structure_based.hmm Data/Processed/neg_${fold}.fasta
    hmmsearch -Z 1000 --max --tblout Data/Processed/seq_pos_${fold}.out \
              Data/Processed/sequence_based.hmm Data/Processed/pos_${fold}.fasta
    hmmsearch -Z 1000 --max --tblout Data/Processed/seq_neg_${fold}.out \
              Data/Processed/sequence_based.hmm Data/Processed/neg_${fold}.fasta
done
```

### Step 5 — Evaluate performance

```bash
# Build .class files (format: ID label E-value; undetected sequences get E-value = 10)
python3 Scripts/build_class_files.py

# Sweep E-value thresholds 10⁻¹ → 10⁻¹²
for model in struct seq; do
    for fold in 1 2; do
        for i in $(seq 1 12); do
            python3 Scripts/performance.py \
                Data/Processed/${model}_set_${fold}.class 1e-$i
        done
    done
done
```

### Step 6 — Generate figures

```bash
jupyter notebook Scripts/MCC_vs_Evalue_Plot.ipynb
jupyter notebook Scripts/Confusion_Matrix_Table.ipynb
jupyter notebook Scripts/Evalue_Distribution_Plot.ipynb
```

---

## 📂 Data Availability

Large files are excluded from this repository due to size constraints. They can be retrieved as follows:

| File                       | Source     | Query                                                                       |
|----------------------------|------------|-----------------------------------------------------------------------------|
| `human_kunitz.fasta`       | UniProtKB  | `(database:pfam PF00014) AND reviewed:true AND organism_id:9606`            |
| `nothuman_kunitz.fasta`    | UniProtKB  | `(database:pfam PF00014) AND reviewed:true NOT organism_id:9606`            |
| `swissprot_notkunitz.fasta`| UniProtKB  | `reviewed:true NOT (database:pfam PF00014)`                                 |
| PDB structures             | RCSB batch | IDs in `Data/Raw/pdb_ids.txt`                                               |

The trained HMMs (`structure_based.hmm`, `sequence_based.hmm`) and both alignments are included in the repository and can be used directly for Kunitz domain detection without re-running the full pipeline.


