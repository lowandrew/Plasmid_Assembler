#!/usr/bin/env python
import multiprocessing
import argparse
import glob
import os
from Bio import SeqIO
from biotools import kmc
from biotools import mash
from biotools import bbtools
from biotools.accessoryfunctions import file_len
from accessoryFunctions import accessoryFunctions
# My previous PlasmidExtractor code was fairly crappy (and not testable!) so it's getting rewritten here.

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

    if reverse_reads:  # If reverse reads file is provided, do things in paired end mode.
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
                                          reverse_out=os.path.join(output_dir, 'plasmid_reads_R1.fastq.gz'),
                                          rskip=6,
                                          threads=threads)
        else:  # If we have a computer with lots of memory, don't use reduced mode.
            out, err = bbtools.bbduk_bait(reference=plasmid_db,
                                          forward_in=os.path.join(output_dir, 'reads_trimmed_R1.fastq.gz'),
                                          reverse_in=os.path.join(output_dir, 'reads_trimmed_R2.fastq.gz'),
                                          forward_out=os.path.join(output_dir, 'plasmid_reads_R1.fastq.gz'),
                                          reverse_out=os.path.join(output_dir, 'plasmid_reads_R1.fastq.gz'),
                                          threads=threads)
        if logfile:  # Write out and err from baiting.
            accessoryFunctions.write_to_logfile(out, err, logfile)
    else:  # Only forward reads means do things in single end mode.
        # First trim.
        out, err = bbtools.bbduk_trim(forward_in=forward_reads,
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
    results = mash.read_mash_screen(screen_result=os.path.join(output_dir, 'screent_results.tsv'))
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
        out, err = kmc.kmc(forward_in=os.path.join(fasta_dir, potential_plasmid_list),
                           database_name=os.path.join(output_dir, plasmid),
                           tmpdir=os.path.join(output_dir, 'tmp'),
                           fm='',
                           t=threads)
        if logfile:
            accessoryFunctions.write_to_logfile(out, err, logfile)


def find_plasmid_kmer_scores(reads_kmerized, kmc_database_dir, output_dir,
                             threads=1):
    score_dict = dict()
    kmerized_plasmids = glob.glob(os.path.join(kmc_database_dir, '*.kmc_pre'))
    read_list = [reads_kmerized] * len(kmerized_plasmids)
    pool = multiprocessing.Pool(processes=threads)
    results = pool.starmap(find_score, read_list, kmerized_plasmids)
    pool.close()
    pool.join()


def find_score(read_db, plasmid_db):  # THIS NEEDS TO BE FINSIHED.
    plasmid = plasmid_db.replace('.kmc_pre', '')
    kmc.intersect(read_db, plasmid_db, plasmid + '_intersect')
    kmc.dump(plasmid, plasmid + '_dump')
    kmc.dump(plasmid_db + '_intersect', plasmid + '_intersect_dump')
    num_plasmid_kmers = file_len(plasmid + '_dump')
    num_intersect_kmers = file_len(plasmid + '_intersect_dump')
    score = float(num_intersect_kmers) / float(num_plasmid_kmers)
    return score


if __name__ == '__main__':
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
    args = parser.parse_args()
