# save as fix_class_files.py
from Bio import SeqIO

def parse_tblout_hits(filepath):
    """Returns dict of protein_id -> evalue for sequences that had hits"""
    hits = {}
    with open(filepath) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            cols = line.split()
            hits[cols[0]] = cols[4]
    return hits

def build_class(pos_fasta, neg_fasta, pos_out, neg_out, outfile):
    pos_hits = parse_tblout_hits(pos_out)
    neg_hits = parse_tblout_hits(neg_out)
    
    lines = []
    
    # All positives (label=1)
    for r in SeqIO.parse(pos_fasta, 'fasta'):
        pid = r.id
        evalue = pos_hits.get(pid, '10')  # no hit → E-value = 10
        lines.append(f"{pid} 1 {evalue}")
    
    # All negatives (label=0)
    for r in SeqIO.parse(neg_fasta, 'fasta'):
        pid = r.id
        evalue = neg_hits.get(pid, '10')  # no hit → E-value = 10
        lines.append(f"{pid} 0 {evalue}")
    
    with open(outfile, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    print(f"Written: {outfile}  ({len(lines)} total entries)")

build_class(
    'Data/Processed/pos_1.fasta', 'Data/Processed/neg_1.fasta',
    'Data/Processed/struct_pos_1.out', 'Data/Processed/struct_neg_1.out',
    'Data/Processed/struct_set_1.class'
)
build_class(
    'Data/Processed/pos_2.fasta', 'Data/Processed/neg_2.fasta',
    'Data/Processed/struct_pos_2.out', 'Data/Processed/struct_neg_2.out',
    'Data/Processed/struct_set_2.class'
)
build_class(
    'Data/Processed/pos_1.fasta', 'Data/Processed/neg_1.fasta',
    'Data/Processed/seq_pos_1.out', 'Data/Processed/seq_neg_1.out',
    'Data/Processed/seq_set_1.class'
)
build_class(
    'Data/Processed/pos_2.fasta', 'Data/Processed/neg_2.fasta',
    'Data/Processed/seq_pos_2.out', 'Data/Processed/seq_neg_2.out',
    'Data/Processed/seq_set_2.class'
)
