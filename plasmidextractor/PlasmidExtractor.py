#!/usr/bin/env python
import multiprocessing
import subprocess
import argparse
import shutil
import glob
import time
import os
from Bio import SeqIO
from biotools import kmc
from scipy import cluster
from biotools import mash
from biotools import bbtools
from biotools.accessoryfunctions import file_len
from accessoryFunctions import accessoryFunctions
from accessoryFunctions.accessoryFunctions import printtime

"""
Methods that we'll need to include in this version of the code:
MAKE EVERYTHING WORK ON BOTH PAIRED AND UNPAIRED READS...

Plasmid Read Baiting/Quality Trimming
Mash Screening For Potential Plasmids
Kmerization of Potential Plasmids
Kmerization of Reference
Plasmid Scoring based on kmer intersection
Plasmid Filtering for Similarity
Detection of Incompatibility/Virulence/AMR genes
"""


def bait_and_trim(forward_reads, output_dir, plasmid_db, logfile=None, reverse_reads=None, low_memory=False,
                  threads=1):
    # If for some reason the output dir hasn't been created, do so now.
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    if reverse_reads is not None:  # If reverse reads file is provided, do things in paired end mode.
        # First trim.
        out, err = bbtools.bbduk_trim(forward_in=forward_reads,
                                      reverse_in=reverse_reads,
                                      forward_out=os.path.join(output_dir, 'reads_trimmed_R1.fastq.gz'),
                                      reverse_out=os.path.join(output_dir, 'reads_trimmed_R2.fastq.gz'),
                                      threads=threads)
        if logfile:  # Write out and err from trimming.
            accessoryFunctions.write_to_logfile(out, err, logfile)
        # Now bait out plasmid reads.
        if low_memory:  # Use less memory if that was specified.
            out, err = bbtools.bbduk_bait(reference=plasmid_db,
                                          forward_in=os.path.join(output_dir, 'reads_trimmed_R1.fastq.gz'),
                                          reverse_in=os.path.join(output_dir, 'reads_trimmed_R2.fastq.gz'),
                                          forward_out=os.path.join(output_dir, 'plasmid_reads_R1.fastq.gz'),
                                          reverse_out=os.path.join(output_dir, 'plasmid_reads_R2.fastq.gz'),
                                          rskip=6,
                                          threads=threads)
        else:  # If we have a computer with lots of memory, don't use reduced mode.
            out, err = bbtools.bbduk_bait(reference=plasmid_db,
                                          forward_in=os.path.join(output_dir, 'reads_trimmed_R1.fastq.gz'),
                                          reverse_in=os.path.join(output_dir, 'reads_trimmed_R2.fastq.gz'),
                                          forward_out=os.path.join(output_dir, 'plasmid_reads_R1.fastq.gz'),
                                          reverse_out=os.path.join(output_dir, 'plasmid_reads_R2.fastq.gz'),
                                          threads=threads)
        if logfile:  # Write out and err from baiting.
            accessoryFunctions.write_to_logfile(out, err, logfile)
    else:  # Only forward reads means do things in single end mode.
        # First trim.
        out, err = bbtools.bbduk_trim(forward_in=forward_reads,
                                      reverse_in='NA',
                                      forward_out=os.path.join(output_dir, 'reads_trimmed_R1.fastq.gz'),
                                      threads=threads)
        if logfile:  # Write out and err from trimming.
            accessoryFunctions.write_to_logfile(out, err, logfile)
        # Now bait out plasmid reads.
        if low_memory:  # Use less memory if that was specified.
            out, err = bbtools.bbduk_bait(reference=plasmid_db,
                                          forward_in=os.path.join(output_dir, 'reads_trimmed_R1.fastq.gz'),
                                          forward_out=os.path.join(output_dir, 'plasmid_reads_R1.fastq.gz'),
                                          rskip=6,
                                          threads=threads)
        else:  # If we have a computer with lots of memory, don't use reduced mode.
            out, err = bbtools.bbduk_bait(reference=plasmid_db,
                                          forward_in=os.path.join(output_dir, 'reads_trimmed_R1.fastq.gz'),
                                          forward_out=os.path.join(output_dir, 'plasmid_reads_R1.fastq.gz'),
                                          threads=threads)
        if logfile:  # Write out and err from baiting.
            accessoryFunctions.write_to_logfile(out, err, logfile)


def mash_for_potential_plasmids(plasmid_db, forward_reads, output_dir, reverse_reads=None, threads=1, logfile=None,
                                identity_cutoff=0.95):
    """
    Uses mash to find a list of potential plasmids in a set of forward (and optionally reverse) reads.
    :param plasmid_db: Path to a multi-Fasta-formatted file that has plasmid sequences of interest.
    :param forward_reads: Path to forward reads.
    :param output_dir: Path to output directory where mash sketch/screen result file will be stored.
    :param reverse_reads: Path to reverse reads. If not specified, things will work in unpaired mode.
    :param threads: Number of threads to run mash analyses on.
    :param logfile: Path to logfile you want to use.
    :param identity_cutoff: Mash screen identity cutoff. Values lower than this won't be reported.
    :return: potential_plasmids: A list where each entry is a putatively present plasmid, identified by
    the fasta header.
    """
    potential_plasmids = list()
    # Make sure the output dir specified gets created if it doesn't exist.
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    # Make a sketch of the plasmid db.
    out, err = mash.sketch(plasmid_db,
                           output_sketch=os.path.join(output_dir, 'plasmid_sketch.msh'),
                           threads=threads,
                           i='')
    if logfile:
        accessoryFunctions.write_to_logfile(out, err, logfile)

    # Now it's time to use mash screen to try to figure out what plasmids might be present in our sample.
    if reverse_reads:  # As usual, do things slightly differently for paired vs unpaired reads.
        out, err = mash.screen(os.path.join(output_dir, 'plasmid_sketch.msh'),
                               forward_reads, reverse_reads,
                               output_file=os.path.join(output_dir, 'screen_results.tsv'),
                               threads=threads,
                               i=identity_cutoff)
        if logfile:
            accessoryFunctions.write_to_logfile(out, err, logfile)
    else:  # Unpaired read mode.
        out, err = mash.screen(os.path.join(output_dir, 'plasmid_sketch.msh'),
                               forward_reads,
                               output_file=os.path.join(output_dir, 'screen_results.tsv'),
                               threads=threads,
                               i=identity_cutoff)
        if logfile:
            accessoryFunctions.write_to_logfile(out, err, logfile)

    # Now need to read through the list of potential plasmids generated by the mash screen.
    results = mash.read_mash_screen(screen_result=os.path.join(output_dir, 'screen_results.tsv'))
    for item in results:
        potential_plasmids.append(item.query_id)

    return potential_plasmids


def create_individual_fastas(plasmid_db, potential_plasmid_list, output_dir):
    """
    Creates individual FASTAs from a multifasta, writes them in output_dir
    :param plasmid_db: Path to multifasta file.
    :param potential_plasmid_list: List of identifiers from multifasta you want to write individual FASTAs
    for.
    :param output_dir: Directory where you want to put the individual FASTAs.
    """
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    for record in SeqIO.parse(plasmid_db, 'fasta'):
        if record.id in potential_plasmid_list:
            SeqIO.write(record, os.path.join(output_dir, record.id), 'fasta')


def kmerize_individual_fastas(potential_plasmid_list, fasta_dir, output_dir, threads=1, logfile=None):
    """
    Creates a KMC database for a list of potential plasmids that have FASTA-formatted sequences in fasta_dir.
    KMC databases are placed in output_dir.
    :param potential_plasmid_list: List of potential plasmids.
    :param fasta_dir: Directory where FASTA files for each potential plasmid are located.
    :param output_dir: Directory to store KMC Databases in. Created if it doesn't exist.
    :param logfile: File to write output to.
    :param threads: Number of threads to run KMC with.
    """
    if not os.path.isdir(output_dir):  # Make output dir if necessary.
        os.makedirs(output_dir)

    for plasmid in potential_plasmid_list:  # Call KMC in FASTA mode on each individual FASTA.
        out, err = kmc.kmc(forward_in=os.path.join(fasta_dir, plasmid),
                           database_name=os.path.join(output_dir, plasmid),
                           tmpdir=os.path.join(output_dir, 'tmp'),
                           fm='',
                           t=threads)
        if logfile:
            accessoryFunctions.write_to_logfile(out, err, logfile)


def find_plasmid_kmer_scores(reads_kmerized, kmc_database_dir, output_dir,
                             threads=1, cutoff=0.95):
    """
    Computes kmer overlaps (how many kmers are in both plasmid and reads) for a set of reads and however many
    kmer databases are in the kmc_database dir folder.
    :param reads_kmerized: Name of kmerized read database.
    :param kmc_database_dir: Folder containing kmc database dirs.
    :param output_dir: Where to store intermediate files.
    :param threads: Number of threads to use.
    :param cutoff: Cutoff for finding
    :return:
    """
    present_plasmids = dict()
    kmerized_plasmids = glob.glob(os.path.join(kmc_database_dir, '*.kmc_pre'))
    if os.path.isfile(os.path.join(output_dir, 'plasmid_reads_R2.fastq.gz')):
        kmc.kmc(forward_in=os.path.join(output_dir, 'plasmid_reads_R1.fastq.gz'),
                reverse_in=os.path.join(output_dir, 'plasmid_reads_R2.fastq.gz'),
                tmpdir=os.path.join(output_dir, 'tmp'),
                database_name=os.path.join(output_dir, 'read_kmers'))
    else:
        kmc.kmc(forward_in=os.path.join(output_dir, 'plasmid_reads_R1.fastq.gz'),
                tmpdir=os.path.join(output_dir, 'tmp'),
                database_name=os.path.join(output_dir, 'read_kmers'))
    read_list = [reads_kmerized] * len(kmerized_plasmids)
    pool = multiprocessing.Pool(processes=threads)
    results = pool.starmap(find_score, zip(read_list, kmerized_plasmids))
    pool.close()
    pool.join()
    for result in results:
        if result[1] > cutoff:
            present_plasmids[result[0]] = result[1]
    return present_plasmids


def find_score(read_db, plasmid_db):
    """
    Given a database of reads created by KMC and a database from a plasmid FASTA created by KMC,
    will find how many plasmid kmers are found in the read dataset.
    :param read_db: KMC database of read kmers.
    :param plasmid_db: KMC database of plasmid kmers.
    :return: List with the name of the plasmid database at index 0, and score for that plasmid at index 1.
    """
    plasmid = plasmid_db.replace('.kmc_pre', '')
    kmc.intersect(read_db, plasmid, plasmid + '_intersect')
    kmc.dump(plasmid, plasmid + '_dump')
    kmc.dump(plasmid + '_intersect', plasmid + '_intersect_dump')
    num_plasmid_kmers = file_len(plasmid + '_dump')
    num_intersect_kmers = file_len(plasmid + '_intersect_dump')
    score = float(num_intersect_kmers) / float(num_plasmid_kmers)
    return [plasmid, score]


def filter_similar_plasmids(plasmid_scores, output_dir):
    """
    Fairly frequently the plasmids recovered end up being close to identical.
    This method sorts through the plasmids to find which plasmids are very similar and picks the best one.
    :param plasmid_scores: Dictionary generated by find_plasmid_kmer_scores.
    :param output_dir: Directory to put temporary files.
    :return: List of plasmids.
    """
    # If only one plasmid, just return that.
    if len(plasmid_scores) == 1:
        return list(plasmid_scores.keys())
    # Otherwise, we create a distance matrix using mash.
    mash_results = list()
    i = 0
    for query_plasmid in plasmid_scores:
        mash_results.append(list())
        for reference_plasmid in plasmid_scores:
            mash.dist(query_plasmid, reference_plasmid, output_file=os.path.join(output_dir, 'distances.tab'))
            result = mash.read_mash_output(os.path.join(output_dir, 'distances.tab'))
            mash_results[i].append(result)
        i += 1
    matrix = list()
    iteration = 1
    for result in mash_results:
        j = 1
        for item in result:
            if j > iteration:
                matrix.append(item[0].distance)
            j += 1
        iteration += 1

    # Once the distance matrix has been made, feed it into SciPy to do clustering.
    z = cluster.hierarchy.linkage(matrix, method='average')
    clustering = cluster.hierarchy.fcluster(z, 0.05, criterion='distance')
    num_clusters = max(clustering)
    clusters = list()
    # Create our clusters.
    plasmids_to_use = list()
    for i in range(num_clusters):
        clusters.append(list())
    plasmid_names = list(plasmid_scores.keys())
    for i in range(len(clustering)):
        clusters[clustering[i] - 1].append(plasmid_names[i])
    # Iterate through clusters, and use the highest scoring plasmid from each cluster for further analysis.
    for group in clusters:
        max_score = 0
        best_hit = ''
        for strain in group:
            if plasmid_scores[strain] > max_score:
                best_hit = strain
                max_score = plasmid_scores[strain]
        plasmids_to_use.append(best_hit)
    return plasmids_to_use


def generate_consensus(forward_reads, reference_fasta, output_fasta, logfile=None,
                       cleanup=True, output_base='out', threads=4, reverse_reads=None, gap_fill=False):
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
    outlog = logfile + '_out.txt'
    errlog = logfile + '_err.txt'
    with open(outlog, 'a+') as outlog_handle:
        with open(errlog, 'a+') as errlog_handle:
            # Step 1: Index fasta file
            cmd = 'samtools faidx {}'.format(reference_fasta)
            subprocess.call(cmd, shell=True, stderr=errlog_handle, stdout=outlog_handle)
            # Step 2: Run bbmap to generate sam/bamfile.
            if reverse_reads:
                cmd = 'bbmap.sh ref={} in={} in2={} out={} nodisk overwrite threads={}'.format(reference_fasta, forward_reads,
                                                                                               reverse_reads, output_base + '.bam',
                                                                                               str(threads))
                subprocess.call(cmd, shell=True, stderr=errlog_handle, stdout=outlog_handle)
            else:
                cmd = 'bbmap.sh ref={} in={} out={} nodisk overwrite threads={}'.format(reference_fasta, forward_reads,
                                                                                        output_base + '.bam',
                                                                                        str(threads))
                subprocess.call(cmd, shell=True, stderr=errlog_handle, stdout=outlog_handle)
            # Step 3: Sort the bam file.
            cmd = 'samtools sort {}.bam -o {}_sorted.bam'.format(output_base, output_base)
            subprocess.call(cmd, shell=True, stderr=errlog_handle, stdout=outlog_handle)
            # Step 3.1: Use bedtools + some shell magic to find regions with zero coverage that we'd like to be Ns.
            # Ideally, change this from awk at some point to make it more generalizable.
            cmd = 'bedtools genomecov -ibam {}_sorted.bam -bga | awk \'$4 == 0\' > {}.bed'.format(output_base, output_base)
            subprocess.call(cmd, shell=True, stderr=errlog_handle, stdout=outlog_handle)
            # Step 3.2: Fancy bcftools piping to generate vcf file.
            cmd = 'bcftools mpileup --threads {} -Ou -f {} {}_sorted.bam | bcftools call  --ploidy ' \
                  '1 --threads {} -mv -Oz -o {}.vcf.gz'\
                  ''.format(str(threads), reference_fasta, output_base, str(threads), output_base)
            subprocess.call(cmd, shell=True, stderr=errlog_handle, stdout=outlog_handle)
            # Step 4: Index vcf file.
            cmd = 'tabix {}.vcf.gz'.format(output_base)
            subprocess.call(cmd, shell=True, stderr=errlog_handle, stdout=outlog_handle)
            # Step 5: Generate consensus fasta from vcf file.
            cmd = 'cat {} | bcftools consensus {}.vcf.gz > {}'.format(reference_fasta, output_base, 'tmp.fasta')
            subprocess.call(cmd, shell=True, stderr=errlog_handle, stdout=outlog_handle)
            # Step 6: Mask regions that don't have any coverage using bedtools.
            cmd = 'bedtools maskfasta -fi tmp.fasta -bed {}.bed -fo {}'.format(output_base, output_fasta)
            subprocess.call(cmd, shell=True, stderr=errlog_handle, stdout=outlog_handle)
            if gap_fill:
                cmd = 'samtools index {}_sorted.bam'.format(output_base)
                subprocess.call(cmd, shell=True, stderr=errlog_handle, stdout=outlog_handle)
                cmd = 'pilon --genome {fasta} --bam {bam} ' \
                      '--fix "all" --output {pilon_fasta}'.format(fasta=output_fasta,
                                                                  bam=output_base + '_sorted.bam',
                                                                  pilon_fasta=output_fasta.replace('.fasta', '_pilon'))
                subprocess.call(cmd, shell=True, stderr=errlog_handle, stdout=outlog_handle)

    if cleanup:
        to_delete = [output_base + '.bam', reference_fasta + '.fai', output_base + '_sorted.bam',
                     output_base + '.vcf.gz', output_base + '.vcf.gz.tbi', output_base + '.bed', 'tmp.fasta']
        for item in to_delete:
            try:
                os.remove(item)
            except FileNotFoundError:
                pass


def do_plasmid_typing(databases_folder, sample_directories):
    """
    Runs plasmid typing for AMR/virulence/incompatibility genes on completed PlasmidExtractor FASTAs.
    :param databases_folder: folder where databases for each kind of typing are stored.
    :param sample_directories: List of directories where PlasmidExtractor output was created.
    """
    for directory in sample_directories:
        # Only try to run GeneSeekr if at least one plasmid was found.
        if len(glob.glob(os.path.join(directory, '*fasta'))) > 0:
            # AMR, then incompatibility, then virulence. Don't need to see out or err for these ones.
            with open(os.devnull, 'w') as out:
                cmd = 'GeneSeekr.py -s {sample_dir} -r {sample_dir} -t {amr_targets} ' \
                      '-v --report_name resistance.csv'.format(sample_dir=directory,
                                                               amr_targets=os.path.join(databases_folder, 'resistance_db'))
                subprocess.call(cmd, shell=True, stderr=out, stdout=out)
                cmd = 'GeneSeekr.py -s {sample_dir} -r {sample_dir} -t {inc_targets} ' \
                      '-v --report_name incompatibility.csv'.format(sample_dir=directory,
                                                                    inc_targets=os.path.join(databases_folder, 'incompatibility'))
                subprocess.call(cmd, shell=True, stderr=out, stdout=out)
                cmd = 'GeneSeekr.py -s {sample_dir} -r {sample_dir} -t {vir_targets} ' \
                      '-v --report_name virulence.csv'.format(sample_dir=directory,
                                                              vir_targets=os.path.join(databases_folder, 'virulence_db'))
                subprocess.call(cmd, shell=True, stderr=out, stdout=out)


def create_summary_report(sample_directories, output_dir, report_name):
    """
    Looks into directories where PlasmidExtractor has run and concatenates reports into one report in output dir.
    :param sample_directories: List of directories where PlasmidExtractor output was created.
    :param output_dir: Directory to put output.
    :param report_name: The report to look for in each sample_directory, including the .csv
    """
    if accessoryFunctions.dependency_check('pilon') is True:
        gap_fill = True
    else:
        gap_fill = False

    with open(os.path.join(output_dir, report_name), 'w') as outfile:
        outfile.write('Strain,Gene,PercentIdentity,PercentCovered,Contig,Location,Sequence\n')
        for directory in sample_directories:
            if os.path.isfile(os.path.join(directory, report_name)):
                with open(os.path.join(directory, report_name)) as infile:
                    lines = infile.readlines()
                    for i in range(1, len(lines)):
                        if gap_fill is True:
                            if 'pilon' in lines[i]:
                                # Replace the first item with whatever the sample name is.
                                to_write = lines[i].split(',')[1:]
                                if len(to_write) > 0:
                                    to_write = ','.join(to_write)
                                    outfile.write(os.path.split(directory)[-1] + ',' + to_write)
                        else:
                            to_write = lines[i].split(',')[1:]
                            if len(to_write) > 0:
                                to_write = ','.join(to_write)
                                outfile.write(os.path.split(directory)[-1] + ',' + to_write)


if __name__ == '__main__':
    start = time.time()
    # Before starting anything, do a check for external dependencies.
    dependencies = ['bbduk.sh', 'bbmap.sh', 'samtools', 'bedtools', 'bcftools',
                    'kmc', 'mash']
    missing_dependencies = list()
    for dependency in dependencies:
        if accessoryFunctions.dependency_check(dependency) is False:
            missing_dependencies.append(dependency)
    if missing_dependencies:
        print('ERROR: The following dependencies are not available on your $PATH: {}\nPlease install them'
              ' and try rerunning PlasmidExtractor.'.format(missing_dependencies))
        quit()

    # Once we've (hopefully) passed the dependency check, do argument parsing and get things going.
    cpu_count = multiprocessing.cpu_count()  # Default threads to use later will use this.
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_directory',
                        required=True,
                        type=str,
                        help='Path to input directory containing your FASTQ files.')
    parser.add_argument('-o', '--output_directory',
                        required=True,
                        type=str,
                        help='Path to output directory where results will be stored. Created if it does not'
                             ' already exist.')
    parser.add_argument('-p', '--plasmid_database',
                        required=True,
                        type=str,
                        help='Path to Plasmid Database.')
    parser.add_argument('-d', '--databases',
                        required=True,
                        type=str,
                        help='Path to folder containing databases for BLAST typing.')
    parser.add_argument('-fid', '--forward_id',
                        type=str,
                        default='_R1',
                        help='Identifier for forward reads when reads are paired.')
    parser.add_argument('-rid', '--reverse_id',
                        type=str,
                        default='_R2',
                        help='Identifier for reverse reads when reads are paired.')
    parser.add_argument('-t', '--threads',
                        type=int,
                        default=cpu_count,
                        help='Number of threads to run analysis on. Default is number of cores on machine.')
    parser.add_argument('-c', '--cutoff',
                        type=float,
                        default=0.99,
                        help='Similarity cutoff for plasmids found in the RefSeq database.')
    parser.add_argument('-l', '--low_memory',
                        default=False,
                        action='store_true',
                        help='When enabled, will use substantially less memory for baiting out plasmid reads.')
    parser.add_argument('-nc', '--no_consensus',
                        default=False,
                        action='store_true',
                        help='When activated, will only find plasmids, not generate sequences and type them.')
    args = parser.parse_args()

    if accessoryFunctions.dependency_check('pilon') is True:
        gap_fill = True
    else:
        gap_fill = False
    printtime('Welcome to PlasmidExtractor! Beginning workflow...', start)
    # Begin by finding both our paired reads and unpaired reads - these have to be treated slightly differently.
    paired_reads = accessoryFunctions.find_paired_reads(args.input_directory,
                                                        forward_id=args.forward_id,
                                                        reverse_id=args.reverse_id)
    unpaired_reads = accessoryFunctions.find_unpaired_reads(args.input_directory,
                                                            forward_id=args.forward_id,
                                                            reverse_id=args.reverse_id)
    # Create our output directory if it doesn't already exist.
    if not os.path.isdir(args.output_directory):
        os.makedirs(args.output_directory)

    # Create a plasmid report that will be appended to as plasmids are found.
    with open(os.path.join(args.output_directory, 'plasmidReport.csv'), 'w') as f:
        f.write('Sample,Plasmid,Score\n')

    sample_directories = list()  # Keep track of sample directories so we know where to look for AMR/Virulence/Inc detection
    # Go through the PlasmidExtractor workflow for paired reads.
    for pair in paired_reads:
        # Get the sample name for the current pair. Use last part of file path, then split on paired read identifier.
        sample_name = os.path.split(pair[0])[-1].split(args.forward_id)[0]
        sample_directories.append(os.path.join(args.output_directory, sample_name))
        log = os.path.join(args.output_directory, sample_name + '.log')
        printtime('Extracting plasmid reads for sample {}...'.format(sample_name), start)
        bait_and_trim(forward_reads=pair[0], reverse_reads=pair[1],
                      output_dir=os.path.join(args.output_directory, sample_name, 'tmp'),
                      plasmid_db=args.plasmid_database, threads=args.threads,
                      low_memory=args.low_memory,
                      logfile=log)
        printtime('Mashing for plasmids...', start)
        plasmids = mash_for_potential_plasmids(forward_reads=os.path.join(args.output_directory, sample_name, 'tmp', 'plasmid_reads_R1.fastq.gz'),
                                               reverse_reads=os.path.join(args.output_directory, sample_name, 'tmp', 'plasmid_reads_R2.fastq.gz'),
                                               output_dir=os.path.join(args.output_directory, sample_name, 'tmp'),
                                               threads=args.threads,
                                               plasmid_db=args.plasmid_database,
                                               identity_cutoff=args.cutoff,
                                               logfile=log)
        printtime('Extracting potential plasmid FASTA files...', start)
        create_individual_fastas(plasmid_db=args.plasmid_database,
                                 potential_plasmid_list=plasmids,
                                 output_dir=os.path.join(args.output_directory, sample_name, 'tmp'))
        printtime('Creating kmer databases...', start)
        kmerize_individual_fastas(potential_plasmid_list=plasmids,
                                  fasta_dir=os.path.join(args.output_directory, sample_name, 'tmp'),
                                  output_dir=os.path.join(args.output_directory, sample_name, 'tmp'),
                                  threads=args.threads,
                                  logfile=log)
        printtime('Scoring plasmids...', start)
        plasmid_scores = find_plasmid_kmer_scores(kmc_database_dir=os.path.join(args.output_directory, sample_name, 'tmp'),
                                                  reads_kmerized=os.path.join(args.output_directory, sample_name, 'tmp', 'read_kmers'),
                                                  output_dir=os.path.join(args.output_directory, sample_name, 'tmp'),
                                                  cutoff=args.cutoff)
        if len(plasmid_scores) == 0:
            printtime('Found zero plasmids...', start)
            continue
        printtime('Filtering out similar plasmids...', start)
        filtered_plasmids = filter_similar_plasmids(plasmid_scores, os.path.join(args.output_directory, sample_name, 'tmp'))
        printtime('Plasmid(s) discovered! Found {} plasmids...'.format(len(filtered_plasmids)), start)
        for plasmid in filtered_plasmids:
            with open(os.path.join(args.output_directory, 'plasmidReport.csv'), 'a+') as f:
                f.write('{sample},{plasmid},{score}\n'.format(sample=sample_name,
                                                              plasmid=os.path.split(plasmid)[-1],
                                                              score=plasmid_scores[plasmid]))
            if args.no_consensus:
                continue

            printtime('Recontructing plasmid similar to {}...'.format(os.path.split(plasmid)[-1]), start)
            generate_consensus(forward_reads=os.path.join(args.output_directory, sample_name, 'tmp', 'plasmid_reads_R1.fastq.gz'),
                               reverse_reads=os.path.join(args.output_directory, sample_name, 'tmp', 'plasmid_reads_R2.fastq.gz'),
                               output_fasta=os.path.join(args.output_directory, sample_name, os.path.split(plasmid)[-1] + '.fasta'),
                               reference_fasta=plasmid,
                               logfile=log,
                               output_base=os.path.join(args.output_directory, sample_name, 'tmp', 'out'),
                               gap_fill=gap_fill)

        shutil.rmtree(os.path.join(args.output_directory, sample_name, 'tmp'))

    # Go through the PlasmidExtractor workflow for unpaired reads.
    for reads in unpaired_reads:
        sample_name = os.path.split(reads)[-1]
        sample_name = os.path.splitext(sample_name)[0]
        sample_directories.append(os.path.join(args.output_directory, sample_name))
        log = os.path.join(args.output_directory, sample_name + '.log')
        printtime('Extracting plasmid reads for sample {}...'.format(sample_name), start)
        bait_and_trim(forward_reads=reads,
                      output_dir=os.path.join(args.output_directory, sample_name, 'tmp'),
                      plasmid_db=args.plasmid_database, threads=args.threads,
                      low_memory=args.low_memory,
                      logfile=log)
        printtime('Mashing for plasmids...', start)
        plasmids = mash_for_potential_plasmids(forward_reads=os.path.join(args.output_directory, sample_name, 'tmp', 'plasmid_reads_R1.fastq.gz'),
                                               output_dir=os.path.join(args.output_directory, sample_name, 'tmp'),
                                               threads=args.threads,
                                               plasmid_db=args.plasmid_database,
                                               identity_cutoff=args.cutoff,
                                               logfile=log)
        printtime('Extracting potential plasmid FASTA files...', start)
        create_individual_fastas(plasmid_db=args.plasmid_database,
                                 potential_plasmid_list=plasmids,
                                 output_dir=os.path.join(args.output_directory, sample_name, 'tmp'))
        printtime('Creating kmer databases...', start)
        kmerize_individual_fastas(potential_plasmid_list=plasmids,
                                  fasta_dir=os.path.join(args.output_directory, sample_name, 'tmp'),
                                  output_dir=os.path.join(args.output_directory, sample_name, 'tmp'),
                                  threads=args.threads,
                                  logfile=log)
        printtime('Scoring plasmids...', start)
        plasmid_scores = find_plasmid_kmer_scores(kmc_database_dir=os.path.join(args.output_directory, sample_name, 'tmp'),
                                                  reads_kmerized=os.path.join(args.output_directory, sample_name, 'tmp', 'read_kmers'),
                                                  output_dir=os.path.join(args.output_directory, sample_name, 'tmp'),
                                                  cutoff=args.cutoff)
        if len(plasmid_scores) == 0:
            printtime('Found zero plasmids...', start)
            continue
        printtime('Filtering out similar plasmids...', start)
        filtered_plasmids = filter_similar_plasmids(plasmid_scores, os.path.join(args.output_directory, sample_name, 'tmp'))
        printtime('Plasmids discovered! Found {} plasmids...'.format(len(filtered_plasmids)), start)
        for plasmid in filtered_plasmids:
            with open(os.path.join(args.output_directory, 'plasmidReport.csv'), 'a+') as f:
                f.write('{sample},{plasmid},{score}\n'.format(sample=sample_name,
                                                              plasmid=os.path.split(plasmid)[-1],
                                                              score=plasmid_scores[plasmid]))

            if args.no_consensus:
                continue
            printtime('Recontructing plasmid similar to {}...'.format(os.path.split(plasmid)[-1]), start)
            generate_consensus(forward_reads=os.path.join(args.output_directory, sample_name, 'tmp', 'plasmid_reads_R1.fastq.gz'),
                               reverse_reads=os.path.join(args.output_directory, sample_name, 'tmp', 'plasmid_reads_R2.fastq.gz'),
                               output_fasta=os.path.join(args.output_directory, sample_name, os.path.split(plasmid)[-1] + '.fasta'),
                               reference_fasta=plasmid,
                               logfile=log,
                               output_base=os.path.join(args.output_directory, sample_name, 'tmp', 'out'),
                               gap_fill=gap_fill)

        shutil.rmtree(os.path.join(args.output_directory, sample_name, 'tmp'))

    # Now FASTA files for each plasmid we've found have been generated - need to get them searched for
    # AMR/virulence/incompatiblity genes.
    if args.no_consensus:
        printtime('PlasmidExtractor workflow complete!', start, '\033[1;32m')
        quit()

    printtime('Performing plasmid typing...', start)
    do_plasmid_typing(databases_folder=args.databases,
                      sample_directories=sample_directories)

    # Now, create summary reports of all typing info that will be placed into the main output directory.
    printtime('Creating summary reports...', start)
    create_summary_report(sample_directories=sample_directories,
                          output_dir=args.output_directory,
                          report_name='resistance.csv')
    create_summary_report(sample_directories=sample_directories,
                          output_dir=args.output_directory,
                          report_name='incompatibility.csv')
    create_summary_report(sample_directories=sample_directories,
                          output_dir=args.output_directory,
                          report_name='virulence.csv')
    printtime('PlasmidExtractor workflow complete!', start, '\033[1;32m')
