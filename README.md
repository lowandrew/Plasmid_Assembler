# Plasmid Assembler

Plasmid Assembler doesn't really assemble plasmids, but instead finds plasmids in raw reads that match
to known plasmids and constructs a consensus sequence for that plasmid.

### Program Dependencies

- samtools >= 1.6
- bcftools >= 1.6
- ncbi-blast >= 2.2.31
- bedtools >= 2.25.0
- BBtools >= 37.23
- mash >= 1.1.1 (2.0 highly recommended.)
- kmc >= 3.0
- python >= 3.5

### Python package requirements

In `requirements.txt`. Use `pip3 install -r requirements.txt` to download and install.

### Databases you'll need to download
- https://figshare.com/s/18de8bdcbba47dbaba41
- From that database, move nucleotideseqs.fa to your current directory, and then
run extract_seqs.py (`python extract_seqs.py`). This will separate the sequences in this file into 
separate fastas in a folder called `plasmid_sequences`
- Once that has completed, run the kmc_plasmids.py script to pre-compute the kmers present in each plasmid
(`python kmc_plasmids.py`). These pre-computed kmers will be stored in a folder called `kmerized_plasmids`

#### Running Plasmid Extractor

`python3 PlasmidExtractor.py -o output_dir -kdb kmerized_plasmids -sdb plasmid_sequences -i read_directory`
Where output_dir is where you want your results stored, kmerized_plasmids is the folder created using kmc_plasmids.py, and 
plasmid_sequences is the folder created using extract_seqs.py.

If no plasmids are found, the output directory will be empty.

If plasmids are found, reads without plasmid sequence in them, as well as a fasta file for each plasmid, 
will be in the output directory.