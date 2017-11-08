# Plasmid Assembler

For more complete instructions, see [documentation](https://lowandrew.github.io/Plasmid_Assembler)

Plasmid Assembler doesn't really assemble plasmids, but instead finds plasmids in raw reads that match
to known plasmids and (optionally) constructs a consensus sequence for that plasmid.

### Program Dependencies

- samtools >= 1.6
- bcftools >= 1.6
- ncbi-blast >= 2.2.31
- bedtools >= 2.25.0
- BBtools >= 37.23
- mash >= 2.0
- kmc >= 3.0
- python >= 3.5

### Python package requirements

In `requirements.txt`. Use `pip3 install -r requirements.txt` to download and install.

### Databases you'll need to download
- https://figshare.com/s/18de8bdcbba47dbaba41
- From that database, you will need the nucleotideseq.fa file.

#### Running Plasmid Extractor

`python3 Extractor.py -o output_dir -sdb plasmid_sequences -i read_directory`

Where `output_dir` is where you want your results stored, `plasmid_sequences` is the path to the nucleotideseq.fa file that was downloaded, and `read_directory` is the path to a folder containing your
FASTQ sequences to be analyzed.


If plasmids are found, reads without plasmid sequence in them, as well as a fasta file for each plasmid, 
will be in subfolders in the output directory. A file called _plasmidReport.csv_ will be created in the output directory.


#### Options

```
usage: Extractor.py [-h] -o OUTPUT_DIR -sdb SEQUENCE_DB [-t THREADS] -i
                    INPUT_DIRECTORY [-k] [-c CUTOFF] [-r REPORT]
                    [-fid FORWARD_ID] [-rid REVERSE_ID] [-nc]

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Output directory where results will be stored.
  -sdb SEQUENCE_DB, --sequence_db SEQUENCE_DB
                        Path to directory containing plasmid sequences.
  -t THREADS, --threads THREADS
                        Number of CPUs to run analysis on. Defaults to number
                        of cores on your machine.
  -i INPUT_DIRECTORY, --input_directory INPUT_DIRECTORY
                        Path to directory containing the paired fastq files to
                        be analyzed.
  -k, --keep_tmpfiles   When specified, will keep the created tmp directory
                        instead of deleting it.
  -c CUTOFF, --cutoff CUTOFF
                        Similarity cutoff for finding plasmids.
  -r REPORT, --report REPORT
                        Name of report to be created by PlasmidExtractor.
  -fid FORWARD_ID, --forward_id FORWARD_ID
                        Identifier for forward reads.
  -rid REVERSE_ID, --reverse_id REVERSE_ID
                        Identifier for forward reads.
  -nc, --no_consensus   When enabled, consenus sequences will not be
                        generated, which saves a fair bit of time. Report of
                        plasmid presence will still be found in
                        plasmidReport.csv

```
