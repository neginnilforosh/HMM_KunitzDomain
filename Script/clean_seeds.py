
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

records = list(SeqIO.parse('Data/Raw/kunitz_domain_seeds.fasta', 'fasta'))

# Keep only sequences between 45-80 aa (valid Kunitz domain length)
clean = [r for r in records if 45 <= len(r.seq) <= 80]
short = [r for r in records if len(r.seq) < 45]

print("Removed sequences:")
for r in short:
    print(f"  {r.id}: {len(r.seq)} aa — {str(r.seq)}")

print(f"\nKept: {len(clean)} sequences")
SeqIO.write(clean, 'Data/Raw/kunitz_domain_seeds_clean.fasta', 'fasta')