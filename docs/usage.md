# Usage With a Pip Install

The basic usage of Plasmid Extractor is fairly simple. As input, you will need to provide:

- a folder containing paired-end reads that you think may have plasmids in them. It is assumed that forward reads contain 'R1' in their name and reverse reads contain 'R2' - see detailed usage for information on how to change this if your reads are named differently. 
- a plasmid database (included with the PlasmidExtractor distribution - see [Installation](installation.md)) 
- resistance/virulence/incompatibility databases (also included with PlasmidExtractor distribution)
- the path to a folder where you would like your output to be place (this folder will be created if it does not exist).

For example, in order to analyze the reads in the directory `/home/user/reads` and place the output into `output` using the default plasmid database, the command would be:
`Extractor.py -i /home/user/reads -o output -sdb databases/plasmid_db.fasta -d databases`

Within the output directory, you will find the following:

- a CSV file called plasmidReport.csv, which shows the plasmids found for each sample, along with a score that shows (approximately) the percent identity to a reference plasmid
- a folder for each sample, which contains a FASTA file for each plasmid found, a FASTA file with the plasmids concatenated
- if there was more than one sample, an image file showing a dendrogram of similarity of the total plasmid content of each sample, as well as a heatmap showing similarity
- a log file for each sample showing the output from each step of the pipeline 
- CSV files describing virulence, AMR, and incompatibility genes found for each sample.


# Usage with Docker 

Coming soon...

# A Few Usage Examples

The input parameters of Plasmid Extractor can be customized to your liking. Here are a few examples of different things that could be done:

- Set identity cutoff to 0.8, to increase sensitivity at the cost of specificity.
`Extractor.py -i /home/user/reads -o output -sdb databases/plasmid_db.fasta -d databases -c 0.8`

- Keep around temporary files created during execution to take a look at afterwards.
`Extractor.py -i /home/user/reads -o output -sdb databases/plasmid_db.fasta -d databases -k`

- Don't generate sequences for plasmids found, for speedy analysis.
`Extractor.py -i /home/user/reads -o output -sdb databases/plasmid_db.fasta -d databases -nc`

- Run in low memory mode. 
`Extractor.py -i /home/user/reads -o output -sdb databases/plasmid_db.fasta -d databases -l`

- Combine low memory mode with no consensus generation for maximum speed.
`Extractor.py -i /home/user/reads -o output -sdb databases/plasmid_db.fasta -d databases -nc -l`


#### Mandatory Arguments

- `-o, --output_dir`: The location for your output folder. Can be anything you want it to. Output folder will be created if it does not exist.
- `-sdb, --sequence_db`: The path to your sequence database. This database should be an uncompressed multi-fasta file. You can create your own custom database
of plasmid sequences, or append any plasmid sequences you want to the supplied database.
- `-i, --input_directory`: The path to your input directory, which contains your paired-end reads to be analyzed. These reads can be uncompressed, or gzip/bzip2 compressed. It is assumed that forward
 reads contain 'R1' and reverse reads contain 'R2'. This assumption can be changed with the `-fid` and `-rid` options.
- `-d, --databases`: The path to databases for virulence/AMR/incompatibility gene detection.

#### Optional Arguments

- `-t, --threads`: Number of CPUs to run your analysis on. By default this is the number of cores on your machine. Every step of the pipeline except for consensus sequence generation parallelizes quite well, so it's recommended that this be left at the default unless you need to run other jobs at the same time.
- `-k, --keep_tmpfiles`: By default, a temporary directory is created in the output folder specified with `-o` for each sample and deleted once analysis is complete. Temporary files such as quality-trimmed reads and plasmid-only reads are kept here. Specifying this option will keep these temporary files instead of deleting them, in case you want to inspect intermediate files more closely. 
- `-c, --cutoff`: The score cutoff for a plasmid to be considered present. Corresponds roughly to the percent identity between the generated plasmid and the reference plasmid. By default, set to 0.98 to only allow plasmids to be found that are quite similar. Lower values will generally return more plasmids, with less specificity. 
- `-r, --report`: Name for report that will be created in `--output_dir`. Defaults to plasmidReport.csv, but can be changed to whatever you like.
- `-fid, --forward_id`: Identifier for forward reads. Defaults to `R1`, but if your forward reads use a different naming scheme like `_1`, specify `-fid _1` to have these recognized as forward reads.
- `-rid, --reverse_id`: Same as `-fid`, but for reverse reads.
- `-nc, --no_consensus`: Finding consensus sequences takes a fair chunk of time. If you want to skip this step and only identify the plasmids present in your sample, add this option. Adding this option means that any post-analysis of your plasmids that would usually take place (AMR detection, etc) will not occur.
- `-l, --low_memory`: Enabling this flag will make the PlasmidExtractor pipeline use substantially less memory (~7GB of RAM at peak usage instead of >20GB). May cause slight drops in sensitivity, but nothing too drastic. 
- `-rp, --remove_plasmid`: This flag will enable creation of plasmid-free reads, so that you can assemble/do downstream analysis on only chromosomal reads. Cannot be used with the `-nc` flag.


