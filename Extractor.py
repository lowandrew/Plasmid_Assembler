#!/usr/bin/env python
import shutil
import multiprocessing
import argparse
import glob
import os
from biotools import bbtools
from biotools import kmc
from biotools import accessoryfunctions

# The original Plasmid Extractor has largely become an unmanageable mess, as the way it works is drastically
# different from what I'd planned. This will attempt to take its place, with a design that actually makes sense.


"""
Workflow:
1) Look through a folder to find paired fastq reads (and unpaired? TBD).
2) For each pair of reads, do quality trimming with bbduk.
3) Extract the plasmid reads from each sample using bbduk.
4) Kmerize the reads with kmc and then find overlapping kmers with kmerized plasmids to determine which plasmids are
present.
5) Filter out plasmids that are essentially duplicates, taking the ones that are the best match out of the set of
similar plasmids.
6) Generate consensus sequences for plasmids that are putatively present in raw reads.
7) Create (optionally?) plasmid-free reads.
8) Characterize plasmids - in particular, look for AMR and Virulence genes present in the consensus FASTAs.
9) Optionally? look at how similar plasmid content is between samples, and potentially create visualization.

Inputs:
Folder of fastq files.

Outputs:
- CSV file showing each plasmid found for each sample, including percent kmer identity, and plasmid characterization.
- Fasta files for each plasmid found.
- Some sort of plasmid content similarity matrix.

To do this:
Have object for each sample that tracks plasmids present, scores, amr genes present, and all that good stuff.
"""


class PlasmidExtractor:
    def __init__(self, args, reads):
        self.forward_reads = reads[0]
        self.reverse_reads = reads[1]
        if not os.path.isdir(args.output_dir):
            os.makedirs(args.output_dir)
        self.output_base = os.path.join(args.output_dir, os.path.split(self.forward_reads)[-1].split('_')[0])
        self.tmpdir = self.output_base + 'tmp'
        forward_base = os.path.split(self.forward_reads)[-1]
        reverse_base = os.path.split(self.reverse_reads)[-1]
        self.forward_trimmed = os.path.join(self.tmpdir, forward_base.replace('.f', '_trimmed.f'))
        self.reverse_trimmed = os.path.join(self.tmpdir, reverse_base.replace('.f', '_trimmed.f'))
        self.forward_plasmid = os.path.join(self.tmpdir, forward_base.replace('.f', '_plasmid.f'))
        self.reverse_plasmid = os.path.join(self.tmpdir, reverse_base.replace('.f', '_plasmid.f'))
        self.threads = args.threads
        self.keep_tmpfiles = args.keep_tmpfiles
        self.logfile = self.output_base + '.log'
        self.kmer_db = args.kmer_db
        self.sequence_db = args.sequence_db
        self.cutoff = args.cutoff
        self.potential_plasmids = dict()
        self.report = args.report

    def main(self):
        if not os.path.isdir(self.tmpdir):
            os.makedirs(self.tmpdir)
        if not os.path.isdir(self.output_base):
            os.makedirs(self.output_base)
        # Get forward and reverse reads trimmed.
        print('Quality trimming input reads for {sample}...'.format(sample=self.forward_reads))
        bbtools.bbduk_trim(forward_in=self.forward_reads, forward_out=self.forward_trimmed,
                           reverse_in=self.reverse_reads, reverse_out=self.reverse_trimmed, overwrite='t')
        # TODO: Un-hardcode nucleotideseq.fa
        # Extract plasmid reads.
        print('Baiting out plasmid reads for {sample}...'.format(sample=self.forward_reads))
        bbtools.bbduk_bait(reference='nucleotideseq.fa', forward_in=self.forward_trimmed,
                           forward_out=self.forward_plasmid, reverse_in=self.reverse_trimmed,
                           reverse_out=self.reverse_plasmid, overwrite='t')
        # See what plasmids have kmer identity above cutoff. Will make self.potential_plasmids into a dict with keys as
        # the potential plasmids, and values as the kmer identity level.
        print('Searching for potential plasmids...')
        self.find_plasmids()

    def find_plasmids(self):
        read_db = os.path.join(self.tmpdir, 'read_db')
        kmc.kmc(self.forward_plasmid, read_db, reverse_in=self.reverse_plasmid)
        plasmid_files = glob.glob(self.kmer_db + '/*.kmc_pre')
        read_list = [read_db] * len(plasmid_files)
        cutoff_list = [self.cutoff] * len(plasmid_files)
        pool = multiprocessing.Pool(processes=self.threads)
        results = pool.starmap(PlasmidExtractor.compare, zip(plasmid_files, read_list, cutoff_list))
        pool.close()
        pool.join()
        for result in results:
            if result[0] != 'NA':
                self.potential_plasmids[result[0]] = result[1]

    def filter_potential_plasmids(self):
        print('Filtering!')

    @staticmethod
    def compare(plasmid, read_db, cutoff):
        plasmid = plasmid.replace('.kmc_pre', '')
        db = os.path.split(plasmid)[-1]
        kmc.intersect(plasmid, read_db, db + '_intersect')
        kmc.dump(plasmid, db + '_dump')
        kmc.dump(db + '_intersect', db + '_intersect_dump')
        plasmid_kmers = accessoryfunctions.file_len(db + '_dump')
        read_kmers = accessoryfunctions.file_len(db + '_intersect_dump')
        percent = float(read_kmers) / float(plasmid_kmers)
        to_remove = glob.glob(db + '*')
        for item in to_remove:
            try:
                os.remove(item)
            except:
                pass
        if percent > cutoff:
            return [plasmid, percent]
        else:
            return ['NA', percent]


def check_dependencies():
    """
    Makes sure all the programs that need to be installed are installed.
    """
    dependencies = ['samtools', 'bcftools', 'bedtools', 'bbmap.sh', 'bbduk.sh', 'mash', 'kmc']
    not_present = list()
    for dep in dependencies:
        is_present = shutil.which(dep)
        if is_present is None:
            not_present.append(dep)
    if len(not_present) > 0:
        raise ModuleNotFoundError('ERROR! Could not find executable(s) for: {}!'.format(not_present))


def find_paired_reads(fastq_directory, forward_id='R1', reverse_id='R2'):
    """
    Looks at a directory to try to find paired fastq files. Should be able to find anything fastq.
    :param fastq_directory: Complete path to directory containing fastq files.
    :param forward_id: Identifier for forward reads. Default R1.
    :param reverse_id: Identifier for reverse reads. Default R2.
    :return: List containing pairs of fastq files, in format [[forward_1, reverse_1], [forward_2, reverse_2]], etc.
    """
    pair_list = list()
    fastq_files = glob.glob(fastq_directory + '/*.f*q*')
    for name in fastq_files:
        if forward_id in name and os.path.isfile(name.replace(forward_id, reverse_id)):
            pair_list.append([name, name.replace(forward_id, reverse_id)])
    return pair_list


if __name__ == '__main__':
    num_cpus = multiprocessing.cpu_count()
    check_dependencies()
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output_dir',
                        type=str,
                        required=True,
                        help='Output directory where results will be stored.')
    parser.add_argument('-kdb', '--kmer_db',
                        type=str,
                        required=True,
                        help='Path to directory containing kmerized plasmids.')
    parser.add_argument('-sdb', '--sequence_db',
                        type=str,
                        required=True,
                        help='Path to directory containing plasmid sequences.')
    parser.add_argument('-t', '--threads',
                        default=num_cpus, type=int,
                        help='Number of CPUs to run analysis on. Defaults to number of cores on your machine.')
    parser.add_argument('-i', '--input_directory',
                        type=str,
                        required=True,
                        help='Path to directory containing the paired fastq files to be analyzed.')
    parser.add_argument('-k', '--keep-tmpfiles',
                        default=False,
                        action='store_true',
                        help='When specified, will keep the created tmp directory instead of deleting it.')
    parser.add_argument('-c', '--cutoff',
                        default=0.98,
                        type=float,
                        help='Similarity cutoff for finding plasmids.')
    parser.add_argument('-r', '--report',
                        default='plasmidReport.csv',
                        type=str,
                        help='Name of report to be created by PlasmidExtractor.')
    args = parser.parse_args()
    paired_reads = find_paired_reads(args.input_directory)
    for pair in paired_reads:
        extractor = PlasmidExtractor(args, pair)
        extractor.main()
