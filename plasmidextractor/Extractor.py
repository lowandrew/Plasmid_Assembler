#!/usr/bin/env python
import multiprocessing
import subprocess
import argparse
import shutil
import time
import glob
import csv
import os
from Bio import SeqIO
from biotools import kmc
from scipy import cluster
from biotools import mash
from biotools import bbtools
from biotools import accessoryfunctions
from accessoryFunctions import accessoryFunctions

"""
Workflow:
1) Look through a folder to find paired fastq reads (and unpaired? TBD).
2) For each pair of reads, do quality trimming with bbduk.
3) Extract the plasmid reads from each sample using bbduk.
3.5) Do a preliminary screen with Mash to find plasmid presence.
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

"""


class PlasmidExtractor:
    def __init__(self, args, reads, start):
        self.start = start
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
        self.forward_no_plasmid = os.path.join(self.output_base, forward_base.replace('.f', '_noplasmid.f'))
        self.reverse_no_plasmid = os.path.join(self.output_base, reverse_base.replace('.f', '_noplasmid.f'))
        self.threads = args.threads
        self.keep_tmpfiles = args.keep_tmpfiles
        self.logfile = self.output_base + '.log'
        self.sequence_db = args.sequence_db
        self.cutoff = args.cutoff
        self.potential_plasmids = dict()
        self.report = os.path.join(args.output_dir, args.report)
        self.consensus_plasmids = list()
        self.no_consensus = args.no_consensus
        self.low_memory = args.low_memory
        self.clean_reads = args.remove_plasmid

    def main(self):
        if not os.path.isdir(self.tmpdir):
            os.makedirs(self.tmpdir)
        if not os.path.isdir(self.output_base):
            os.makedirs(self.output_base)
        # Get forward and reverse reads trimmed.
        accessoryFunctions.printtime('Quality trimming input reads for {sample}...'.format(sample=self.forward_reads),
                                     self.start)
        out, err = bbtools.bbduk_trim(forward_in=self.forward_reads, forward_out=self.forward_trimmed,
                                      reverse_in=self.reverse_reads, reverse_out=self.reverse_trimmed, overwrite='t')

        with open(self.logfile, 'a+') as logfile:
            logfile.write(out + '\n')
            logfile.write(err + '\n')
        # Extract plasmid reads.
        accessoryFunctions.printtime('Baiting out plasmid reads for {sample}...'.format(sample=self.forward_reads),
                                     self.start)
        if self.low_memory:
            out, err = bbtools.bbduk_bait(reference=self.sequence_db, forward_in=self.forward_trimmed,
                                          forward_out=self.forward_plasmid, reverse_in=self.reverse_trimmed,
                                          reverse_out=self.reverse_plasmid, overwrite='t', rskip='6')
        else:
            out, err = bbtools.bbduk_bait(reference=self.sequence_db, forward_in=self.forward_trimmed,
                                          forward_out=self.forward_plasmid, reverse_in=self.reverse_trimmed,
                                          reverse_out=self.reverse_plasmid, overwrite='t')
        with open(self.logfile, 'a+') as logfile:
            logfile.write(out + '\n')
            logfile.write(err + '\n')
        # See what plasmids have kmer identity above cutoff. Will make self.potential_plasmids into a dict with keys as
        # the potential plasmids, and values as the kmer identity level.
        out, err = mash.sketch(self.sequence_db, i='',
                               threads=self.threads,
                               output_sketch=os.path.join(self.tmpdir, 'plasmid.msh'))
        with open(self.logfile, 'a+') as logfile:
            logfile.write(out + '\n')
            logfile.write(err + '\n')
        # Use mash screen to get a quick list of potential plasmids.
        accessoryFunctions.printtime('Searching reads for plasmids...', self.start)
        out, err = mash.screen(os.path.join(self.tmpdir, 'plasmid.msh'),
                               self.forward_plasmid, self.reverse_plasmid,
                               threads=self.threads,
                               i=self.cutoff, output_file=os.path.join(self.tmpdir, 'screen.tab'))
        with open(self.logfile, 'a+') as logfile:
            logfile.write(out + '\n')
            logfile.write(err + '\n')
        results = mash.read_mash_screen(os.path.join(self.tmpdir, 'screen.tab'))
        # Get a list of potential plasmids from mash screen output.
        plasmids_present = list()
        for result in results:
            plasmids_present.append(result.query_id)
        # Go through the multifasta and write each of the potential plasmids to the tmp directory.
        # Also kmerize the plasmids.
        for record in SeqIO.parse(self.sequence_db, 'fasta'):
            if record.id in plasmids_present:
                SeqIO.write(record, os.path.join(self.tmpdir, record.id), 'fasta')
                # TODO: Get this step parallelized for when we have lots of potential plasmids found by MASH
                kmc.kmc(os.path.join(self.tmpdir, record.id), os.path.join(self.tmpdir, record.id), fm='')
        self.find_plasmids()
        # If no potential plasmids were found, we can skip over the rest of this.
        if len(self.potential_plasmids) > 0:
            # Get a list of plasmids we actually want to use - Often, there are a lot of hits to very similar plasmids.
            # We'll take only the best hit.
            plasmids_to_use = self.filter_potential_plasmids()
            accessoryFunctions.printtime('Found {} plasmids...'.format(str(len(plasmids_to_use))), self.start)
            # Append to our plasmid report, with info on strain, plasmids present, and the similarity score.
            with open(self.report, 'a+') as f:
                for plasmid in plasmids_to_use:
                    f.write('{},{},{}\n'.format(self.output_base, plasmid, str(self.potential_plasmids[plasmid])))
            if not args.no_consensus:
                accessoryFunctions.printtime('Generating consensus sequences...', self.start)
                # Generate consensus plasmids, and also keep track of the file locations for those plasmids, as they get
                # used later.
                for plasmid in plasmids_to_use:
                    generate_consensus(self.forward_plasmid, self.reverse_plasmid, plasmid,
                                       self.output_base + '/' + os.path.split(plasmid)[-1] + '_consensus.fasta',
                                       threads=self.threads, logfile=self.logfile)
                    self.consensus_plasmids.append(self.output_base + '/' + os.path.split(plasmid)[-1] + '_consensus.fasta')
                accessoryFunctions.printtime('Generating plasmid-free reads...', self.start)
                # Remove plasmid from reads.
                if self.clean_reads:
                    self.remove_plasmid_from_reads()

        # Remove temporary files, unless use for some reason wanted to keep them.
        if not self.keep_tmpfiles:
            shutil.rmtree(self.tmpdir)

    def find_plasmids(self):
        read_db = os.path.join(self.tmpdir, 'read_db')
        plasmid_files = glob.glob(self.tmpdir + '/*.kmc_pre')
        kmc.kmc(self.forward_plasmid, read_db, reverse_in=self.reverse_plasmid)
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
        # Having only one potential plasmid causes the attempt a clustering to fail, as you can't really dendrogram
        # with only one item. If this is the case, just return the one plasmid.
        if len(self.potential_plasmids) == 1:
            plasmids_to_use = list(self.potential_plasmids.keys())
            return plasmids_to_use
        # Step 1 in the filtering process: Use mash to find distances between all plasmids.
        mash_results = list()
        i = 0
        for query_plasmid in self.potential_plasmids:
            mash_results.append(list())
            for reference_plasmid in self.potential_plasmids:
                mash.dist(query_plasmid, reference_plasmid, output_file=os.path.join(self.tmpdir, 'distances.tab'))
                x = mash.read_mash_output(os.path.join(self.tmpdir, 'distances.tab'))
                mash_results[i].append(x)
            i += 1
        # Now have a list of all mash results, so pairwise distances are known. Now need to transform the data
        # into a matrix that SciPy can use.
        matrix = list()
        iteration = 1
        for result in mash_results:
            j = 1
            for item in result:
                if j > iteration:
                    matrix.append(item[0].distance)
                j += 1
            iteration += 1
        # Perform UPGMA clustering on the matrix created from mash distance values.
        z = cluster.hierarchy.linkage(matrix, method='average')
        # Get a list of clusters at our cutoff. Testing seems to show that the cutoff at 0.1 works well.
        # Will need to do a somewhat more extensive set of testing in the not too distant future.
        clustering = cluster.hierarchy.fcluster(z, 0.05, criterion='distance')
        num_clusters = max(clustering)
        clusters = list()
        # Create our clusters.
        plasmids_to_use = list()
        for i in range(num_clusters):
            clusters.append(list())
        plasmid_names = list(self.potential_plasmids.keys())
        for i in range(len(clustering)):
            clusters[clustering[i] - 1].append(plasmid_names[i])
        # Iterate through clusters, and use the highest scoring plasmid from each cluster for further analysis.
        for group in clusters:
            max_score = 0
            best_hit = ''
            for strain in group:
                if self.potential_plasmids[strain] > max_score:
                    best_hit = strain
                    max_score = self.potential_plasmids[strain]
            plasmids_to_use.append(best_hit)
        return plasmids_to_use

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

    def remove_plasmid_from_reads(self):
        # First, generate a concatenated fasta.
        with open(os.path.join(self.output_base, 'plasmids_concatenated.fasta'), 'w') as outfile:
            for sequence in self.consensus_plasmids:
                with open(sequence) as infile:
                    outfile.write(infile.read())
        # Filter out reads that contain sequence that is in the concatenated fasta file.
        out, err = bbtools.bbduk_filter(reference=os.path.join(self.output_base, 'plasmids_concatenated.fasta'),
                                        forward_in=self.forward_trimmed,
                                        reverse_in=self.reverse_trimmed, forward_out=self.forward_no_plasmid,
                                        reverse_out=self.reverse_no_plasmid, overwrite='t')
        with open(self.logfile, 'a+') as logfile:
            logfile.write(out + '\n')
            logfile.write(err + '\n')


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
        raise ModuleNotFoundError('ERROR! Could not find executable(s) for: {}!'.format(str(not_present)))


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


def generate_consensus(forward_reads, reverse_reads, reference_fasta, output_fasta, logfile='logfile.log',
                       cleanup=True, output_base='out', threads=4):
    """
    Generates a consensus fasta give a set of interleaved reads and a reference sequence.
    :param forward_reads: Forward input reads.
    :param reverse_reads: Reverse input reads.
    :param reference_fasta: Reference you want to map your reads to.
    :param output_fasta: Output file for your consensus fasta.
    :param logfile: Log to write stdout and stderr to.
    :param cleanup: If true, deletes intermediate files.
    :param output_base: Base name for output files.
    :param threads: Number of threads to run analysis with.
    """
    with open(logfile, 'a+') as log:
        # Step 1: Index fasta file
        cmd = 'samtools faidx {}'.format(reference_fasta)
        subprocess.call(cmd, shell=True, stderr=log, stdout=log)
        # Step 2: Run bbmap to generate sam/bamfile.
        cmd = 'bbmap.sh ref={} in={} in2={} out={} nodisk overwrite threads={}'.format(reference_fasta, forward_reads,
                                                                                       reverse_reads, output_base + '.bam',
                                                                                       str(threads))
        subprocess.call(cmd, shell=True, stderr=log, stdout=log)
        # Step 3: Sort the bam file.
        cmd = 'samtools sort {}.bam -o {}_sorted.bam'.format(output_base, output_base)
        subprocess.call(cmd, shell=True, stderr=log, stdout=log)
        # Step 3.1: Use bedtools + some shell magic to find regions with zero coverage that we'd like to be Ns.
        # Ideally, change this from awk at some point to make it more generalizable.
        cmd = 'bedtools genomecov -ibam {}_sorted.bam -bga | awk \'$4 == 0\' > {}.bed'.format(output_base, output_base)
        subprocess.call(cmd, shell=True, stderr=log, stdout=log)
        # Step 3.2: Fancy bcftools piping to generate vcf file.
        cmd = 'bcftools mpileup --threads {} -Ou -f {} {}_sorted.bam | bcftools call --threads {} -mv -Oz -o {}.vcf.gz'\
              ''.format(str(threads), reference_fasta, output_base, str(threads), output_base)
        subprocess.call(cmd, shell=True, stderr=log, stdout=log)
        # Step 4: Index vcf file.
        cmd = 'tabix {}.vcf.gz'.format(output_base)
        subprocess.call(cmd, shell=True, stderr=log, stdout=log)
        # Step 5: Generate consensus fasta from vcf file.
        cmd = 'cat {} | bcftools consensus {}.vcf.gz > {}'.format(reference_fasta, output_base, 'tmp.fasta')
        subprocess.call(cmd, shell=True, stderr=log, stdout=log)
        # Step 6: Mask regions that don't have any coverage using bedtools.
        cmd = 'bedtools maskfasta -fi tmp.fasta -bed {}.bed -fo {}'.format(output_base, output_fasta)
        subprocess.call(cmd, shell=True, stderr=log, stdout=log)

    if cleanup:
        to_delete = [output_base + '.bam', reference_fasta + '.fai', output_base + '_sorted.bam',
                     output_base + '.vcf.gz', output_base + '.vcf.gz.tbi', output_base + '.bed', 'tmp.fasta']
        for item in to_delete:
            try:
                os.remove(item)
            except FileNotFoundError:
                pass


class PlasmidAnalyzer:
    def __init__(self, args, start):
        self.sourmash_matrix = os.path.join(args.output_dir, 'plasmid')
        self.output_dir = args.output_dir
        self.keep_tmpfiles = args.keep_tmpfiles
        self.logfile = self.output_dir + '/sourmash.log'
        self.databases = args.databases
        self.start = start

    def sourmash(self):
        concatenated_fastas = glob.glob(self.output_dir + '/*/plasmids_concatenated.fasta')
        # Don't bother if there's only one sample, because then there's nothing to compare to.
        if len(concatenated_fastas) > 1:
            cmd = 'sourmash compute --force --output {} '.format(os.path.join(self.output_dir, 'computed'))
            for fasta in concatenated_fastas:
                cmd += fasta + ' '
            with open(self.logfile, 'w') as logfile:
                subprocess.call(cmd, shell=True, stdout=logfile, stderr=logfile)
                cmd = 'sourmash compare --output {} {}'.format(os.path.join(self.output_dir, 'compared'),
                                                               os.path.join(self.output_dir, 'computed'))
                subprocess.call(cmd, shell=True, stdout=logfile, stderr=logfile)
                # Can't find any way to specify output location for output file. To be investigated.
                cmd = 'sourmash plot --labels {}'.format(os.path.join(self.output_dir, 'compared'))
                subprocess.call(cmd, shell=True, stdout=logfile, stderr=logfile)

            if not self.keep_tmpfiles:
                os.remove(os.path.join(self.output_dir, 'compared'))
                os.remove(os.path.join(self.output_dir, 'computed'))

            # Move dendrogram/matrix png files to results folder, since sourmash refuses to have the output name
            # as an option.
            # Leave this try/except until I've actually tested that this works.
            try:
                shutil.move('compared.matrix.png', os.path.join(self.output_dir, 'compared.matrix.png'))
                shutil.move('compared.dendro.png', os.path.join(self.output_dir, 'compared.dendro.png'))
            except:
                pass

    def amr_detection(self):
        accessoryFunctions.printtime('Finding resistance genes...', self.start, '\033[1;35m')
        # Get dictionary of what genes correspond to what resistances from the notes.txt in the resistance_db folder
        resistance_dict = dict()
        with open(os.path.join(self.output_dir, 'resistance.csv'), 'w') as f:
            f.write('Sample,Plasmid,Gene,Resistance,Coverage,PercentIdentity\n')
        with open(os.path.join(self.databases, 'resistance_db/notes.txt')) as infile:
            lines = infile.readlines()
        for line in lines:
            x = line.split(':')
            if '#' not in line:
                resistance_dict[x[0]] = x[1]
        # Find a list of all your samples.
        samples = glob.glob(os.path.join(self.output_dir, '*/'))
        for sample in samples:
            # Run GeneSeekr using CGE's resistance database.
            cmd = 'GeneSeekr.py -t {databases}/resistance_db -s {sequence_folder} ' \
                  '-r {sequence_folder} -v'.format(sequence_folder=sample, databases=self.databases)
            with open(os.devnull, 'w') as null:
                subprocess.call(cmd, shell=True, stdout=null, stderr=null)
            try:  # Rename the output from the GeneSeekr default
                os.rename('{sequence_folder}/virulence.csv'.format(sequence_folder=sample),
                          '{sequence_folder}/resistance.csv'.format(sequence_folder=sample))
            except FileNotFoundError:
                pass
            try:  # Write a summary report for all samples to be placed in the output folder.
                with open('{sequence_folder}/resistance.csv'.format(sequence_folder=sample)) as csvfile:
                    reader = csv.DictReader(csvfile)
                    for row in reader:
                        if row['Gene'] is not None:
                            with open(os.path.join(self.output_dir, 'resistance.csv'), 'a+') as f:
                                try:
                                    resistance = resistance_dict[row['Gene']]
                                except KeyError:
                                    resistance = 'Unknown'
                                f.write('{sample},{plasmid},{gene},{resistance},{coverage},'
                                        '{percent}\n'.format(sample=os.path.split(sample)[-2],
                                                             plasmid=row['Contig'],
                                                             gene=row['Gene'],
                                                             resistance=resistance,
                                                             coverage=row['PercentCovered'],
                                                             percent=row['PercentIdentity']))
                to_remove = glob.glob(os.path.join(sample, '*csv'))  # Remove temporary csv files created by GeneSeekr.
                for item in to_remove:
                    if 'resistance' not in item:
                        os.remove(item)
            except FileNotFoundError:
                pass

    def virulence_detection(self):
        accessoryFunctions.printtime('Finding virulence genes...', self.start, '\033[1;35m')  # Probably change this away from purple.
        with open(os.path.join(self.output_dir, 'virulence.csv'), 'w') as f:
            f.write('Sample,Plasmid,Gene,Coverage,PercentIdentity\n')
        samples = glob.glob(os.path.join(self.output_dir, '*/'))  # Find our samples.
        for sample in samples:
            # Run GeneSeekr on CGE's virulence DB.
            cmd = 'GeneSeekr.py -t {databases}/virulence_db -s {sequence_folder} ' \
                  '-r {sequence_folder} -v'.format(sequence_folder=sample, databases=self.databases)
            with open(os.devnull, 'w') as null:
                subprocess.call(cmd, shell=True, stdout=null, stderr=null)
            try:
                # Write results to the summary for all samples
                with open('{sequence_folder}/virulence.csv'.format(sequence_folder=sample)) as csvfile:
                    reader = csv.DictReader(csvfile)
                    for row in reader:
                        if row['Gene'] is not None:
                            with open(os.path.join(self.output_dir, 'virulence.csv'), 'a+') as f:
                                f.write('{sample},{plasmid},{gene},{coverage},{percent}\n'.format(sample=os.path.split(sample)[-2],  # Clean this up, displays too much.
                                                                                                  plasmid=row['Contig'],
                                                                                                  gene=row['Gene'],
                                                                                                  coverage=row['PercentCovered'],
                                                                                                  percent=row['PercentIdentity']))
                to_remove = glob.glob(os.path.join(sample, '*csv'))  # Get rid of GeneSeekr files.
                for item in to_remove:
                    if 'resistance' not in item and 'virulence' not in item and 'incompatibility' not in item:
                        os.remove(item)
            except FileNotFoundError:
                pass

    def inc_detection(self):
        accessoryFunctions.printtime('Finding incompatibility groups...', self.start, '\033[1;35m')
        # Get dictionary of what genes correspond to what resistances from the notes.txt in the resistance_db folder
        with open(os.path.join(self.output_dir, 'incompatibility.csv'), 'w') as f:
            f.write('Sample,Plasmid,Gene,Coverage,PercentIdentity\n')
        samples = glob.glob(os.path.join(self.output_dir, '*/'))
        for sample in samples:
            cmd = 'GeneSeekr.py -t {databases}/incompatibility -s {sequence_folder} ' \
                  '-r {sequence_folder} -v'.format(sequence_folder=sample, databases=self.databases)
            with open(os.devnull, 'w') as null:
                subprocess.call(cmd, shell=True, stdout=null, stderr=null)
            try:
                os.rename('{sequence_folder}/virulence.csv'.format(sequence_folder=sample),
                          '{sequence_folder}/incompatibility.csv'.format(sequence_folder=sample))
            except FileNotFoundError:
                pass
            try:
                with open('{sequence_folder}/incompatibility.csv'.format(sequence_folder=sample)) as csvfile:
                    reader = csv.DictReader(csvfile)
                    for row in reader:
                        if row['Gene'] is not None:
                            with open(os.path.join(self.output_dir, 'incompatibility.csv'), 'a+') as f:
                                f.write('{sample},{plasmid},{gene},{coverage},{percent}\n'.format(sample=os.path.split(sample)[-2],  # Clean this up, displays too much.
                                                                                                  plasmid=row['Contig'],
                                                                                                  gene=row['Gene'],
                                                                                                  coverage=row['PercentCovered'],
                                                                                                  percent=row['PercentIdentity']))
                to_remove = glob.glob(os.path.join(sample, '*csv'))
                for item in to_remove:
                    if 'resistance' not in item and 'incompatibility' not in item:
                        os.remove(item)
            except FileNotFoundError:
                pass

    def main(self):
        self.sourmash()
        self.amr_detection()
        self.inc_detection()
        self.virulence_detection()


if __name__ == '__main__':
    start = time.time()
    num_cpus = multiprocessing.cpu_count()
    check_dependencies()
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output_dir',
                        type=str,
                        required=True,
                        help='Output directory where results will be stored.')
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
    parser.add_argument('-k', '--keep_tmpfiles',
                        default=False,
                        action='store_true',
                        help='When specified, will keep the created tmp directory instead of deleting it.')
    parser.add_argument('-c', '--cutoff',
                        default=0.99,
                        type=float,
                        help='Similarity cutoff for finding plasmids.')
    parser.add_argument('-r', '--report',
                        default='plasmidReport.csv',
                        type=str,
                        help='Name of report to be created by PlasmidExtractor.')
    parser.add_argument('-fid', '--forward_id',
                        default='R1',
                        type=str,
                        help='Identifier for forward reads.')
    parser.add_argument('-rid', '--reverse_id',
                        default='R2',
                        type=str,
                        help='Identifier for forward reads.')
    parser.add_argument('-nc', '--no_consensus',
                        default=False,
                        action='store_true',
                        help='When enabled, consenus sequences will not be generated, which saves a fair bit of time.'
                             ' Report of plasmid presence will still be found in plasmidReport.csv')
    parser.add_argument('-d', '--databases',
                        required=True,
                        type=str,
                        help='Path to resistance/virulence/plasmid typing databases.')
    parser.add_argument('-l', '--low_memory',
                        default=False,
                        action='store_true',
                        help='When enabled, will use substantially less memory (~ 7GB instead of ~24GB). May come at'
                             ' the cost of some sensitivity.')
    parser.add_argument('-rp', '--remove_plasmid',
                        default=False,
                        action='store_true',
                        help='If enabled, will remove plasmid reads from the input reads specified.')
    args = parser.parse_args()
    # Get a logfile set up.
    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir)
    with open(os.path.join(args.output_dir, args.report), 'w') as f:
        f.write('Sample,Plasmid,Score\n')
    paired_reads = find_paired_reads(args.input_directory, forward_id=args.forward_id, reverse_id=args.reverse_id)
    for pair in paired_reads:
        accessoryFunctions.printtime('Beginning plasmid extraction for {}...'.format(pair[0]), start)
        extractor = PlasmidExtractor(args, pair, start)
        try:
            extractor.main()
        except subprocess.CalledProcessError:
            accessoryFunctions.printtime('Error encountered for {}. Skipping...'.format(pair[0]), start, '\033[1;31m')
            pass
        accessoryFunctions.printtime('Finished plasmid extraction for {}...'.format(pair[0]), start)
    # Now do some post-analysis. Sourmash for nice visualization, AMR gene finding, and maybe more?...
    if not args.no_consensus:
        analyzer = PlasmidAnalyzer(args, start)
        analyzer.main()
