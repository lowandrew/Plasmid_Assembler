from biotools import kmc
import os
import glob


if not os.path.isdir('kmerized_plasmids'):
    os.makedirs('kmerized_plasmids')

plasmids = glob.glob('plasmid_sequences/*')
for plasmid in plasmids:
    kmc.kmc(plasmid, plasmid.replace('plasmid_sequences', 'kmerized_plasmids'))

