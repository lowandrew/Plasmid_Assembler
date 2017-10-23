from Bio import SeqIO
import os


if not os.path.isdir('plasmid_sequences'):
    os.makedirs('plasmid_sequences')

for record in SeqIO.parse('nucleotideseq.fa', 'fasta'):
    f = open('plasmid_sequences/' + record.id, 'w')
    f.write('>' + record.id + '\n')
    f.write(str(record.seq) + '\n')
    f.close()
