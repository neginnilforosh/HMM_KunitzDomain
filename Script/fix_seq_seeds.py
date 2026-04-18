
from Bio import SeqIO, AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Read the repo's structural alignment (contains only true Kunitz chains)
aln = AlignIO.read('./HMM_KunitzDomain-main/Data/Processed/pdb_kunitz_nr_clean.ali', 'fasta')

# Strip gaps to get raw sequences
clean_records = []
for record in aln:
    ungapped = str(record.seq).replace('-', '').replace('.', '')
    if len(ungapped) >= 30:  # sanity check
        clean_records.append(SeqRecord(Seq(ungapped), id=record.id, description=''))

SeqIO.write(clean_records, 'kunitz_true_seeds.fasta', 'fasta')
print(f"Extracted {len(clean_records)} true Kunitz domain sequences:")
for r in clean_records:
    print(f"  {r.id}: {len(r.seq)} aa")
