
from Bio import SeqIO

# Parse domain table to get coordinates
domains = {}
with open('Data/Processed/domain_hits.tbl') as f:
    for line in f:
        if line.startswith('#') or not line.strip():
            continue
        cols = line.split()
        seq_id   = cols[0]
        evalue   = float(cols[12])   # domain i-evalue
        env_from = int(cols[19])     # envelope start (1-based)
        env_to   = int(cols[20])     # envelope end
        if evalue > 1e-3:
            continue
        # Keep best (lowest evalue) domain per sequence
        if seq_id not in domains or evalue < domains[seq_id][0]:
            domains[seq_id] = (evalue, env_from, env_to)

print(f'Found domain coordinates for {len(domains)} sequences')

# Extract domain subsequences
seqs = SeqIO.to_dict(SeqIO.parse('Data/Processed/ok_kunitz_clean.fasta', 'fasta'))
extracted = []
skipped = 0
for seq_id, (evalue, start, end) in domains.items():
    if seq_id not in seqs:
        skipped += 1
        continue
    full_seq = seqs[seq_id]
    domain_seq = full_seq[start-1:end]  # convert to 0-based
    domain_seq.id = seq_id
    domain_seq.description = f'Kunitz_domain_{start}-{end}'
    if 40 <= len(domain_seq) <= 100:   # sanity check
        extracted.append(domain_seq)

SeqIO.write(extracted, 'Data/Processed/kunitz_domains_swssprot.fasta', 'fasta')
print(f'Extracted {len(extracted)} domain sequences (skipped {skipped})')

# Check length distribution
lengths = [len(r.seq) for r in extracted]
print(f'Length range: {min(lengths)}-{max(lengths)} aa, mean={sum(lengths)/len(lengths):.1f}')