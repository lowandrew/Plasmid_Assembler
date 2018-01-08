import os
import shutil

"""
Remaining things to test:
find_plasmid_kmer_scores
find_score
filter_similar_plasmids
generate_consensus
"""
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(0, parentdir)

from New_Extractor import *


def test_mash_paired_gzipped():
    mash_for_potential_plasmids(forward_reads='tests/test_fastqs/paired_R1.fastq.gz',
                                reverse_reads='tests/test_fastqs/paired_R2.fastq.gz',
                                plasmid_db='tests/test_fasta/dummy_db.fasta',
                                output_dir='tests/mash',
                                identity_cutoff=-1)
    assert os.path.isfile('tests/mash/screen_results.tsv')
    shutil.rmtree('tests/mash')


def test_mash_unpaired_gzipped():
    mash_for_potential_plasmids(forward_reads='tests/test_fastqs/paired_R1.fastq.gz',
                                plasmid_db='tests/test_fasta/dummy_db.fasta',
                                output_dir='tests/mash',
                                identity_cutoff=-1)
    assert os.path.isfile('tests/mash/screen_results.tsv')
    shutil.rmtree('tests/mash')


def test_mash_paired_uncompressed():
    mash_for_potential_plasmids(forward_reads='tests/test_fastqs/paired_R1.fastq',
                                reverse_reads='tests/test_fastqs/paired_R2.fastq',
                                plasmid_db='tests/test_fasta/dummy_db.fasta',
                                output_dir='tests/mash',
                                identity_cutoff=-1)
    assert os.path.isfile('tests/mash/screen_results.tsv')
    shutil.rmtree('tests/mash')


def test_mash_unpaired_uncompressed():
    mash_for_potential_plasmids(forward_reads='tests/test_fastqs/paired_R1.fastq',
                                plasmid_db='tests/test_fasta/dummy_db.fasta',
                                output_dir='tests/mash',
                                identity_cutoff=-1)
    assert os.path.isfile('tests/mash/screen_results.tsv')
    shutil.rmtree('tests/mash')


def test_bait_and_trim_paired_gzipped():
    bait_and_trim(forward_reads='tests/test_fastqs/paired_R1.fastq.gz',
                  reverse_reads='tests/test_fastqs/paired_R2.fastq.gz',
                  plasmid_db='tests/test_fasta/dummy_db.fasta',
                  output_dir='tests/out')
    assert os.path.isfile('tests/out/plasmid_reads_R1.fastq.gz') and os.path.isfile('tests/out/plasmid_reads_R2.fastq.gz')
    shutil.rmtree('tests/out')


def test_bait_and_trim_paired_uncompressed():
    bait_and_trim(forward_reads='tests/test_fastqs/paired_R1.fastq',
                  reverse_reads='tests/test_fastqs/paired_R2.fastq',
                  plasmid_db='tests/test_fasta/dummy_db.fasta',
                  output_dir='tests/out')
    assert os.path.isfile('tests/out/plasmid_reads_R1.fastq.gz') and os.path.isfile('tests/out/plasmid_reads_R2.fastq.gz')
    shutil.rmtree('tests/out')


def test_bait_and_trim_unpaired_gzipped():
    bait_and_trim(forward_reads='tests/test_fastqs/paired_R1.fastq.gz',
                  plasmid_db='tests/test_fasta/dummy_db.fasta',
                  output_dir='tests/out')
    assert os.path.isfile('tests/out/plasmid_reads_R1.fastq.gz')
    shutil.rmtree('tests/out')


def test_bait_and_trim_unpaired_uncompressed():
    bait_and_trim(forward_reads='tests/test_fastqs/paired_R1.fastq',
                  plasmid_db='tests/test_fasta/dummy_db.fasta',
                  output_dir='tests/out')
    assert os.path.isfile('tests/out/plasmid_reads_R1.fastq.gz')
    shutil.rmtree('tests/out')


def test_bait_and_trim_paired_gzipped_lowmem():
    bait_and_trim(forward_reads='tests/test_fastqs/paired_R1.fastq.gz',
                  reverse_reads='tests/test_fastqs/paired_R2.fastq.gz',
                  plasmid_db='tests/test_fasta/dummy_db.fasta',
                  output_dir='tests/out',
                  low_memory=True)
    assert os.path.isfile('tests/out/plasmid_reads_R1.fastq.gz') and os.path.isfile('tests/out/plasmid_reads_R2.fastq.gz')
    shutil.rmtree('tests/out')


def test_bait_and_trim_unpaired_gzipped_lowmem():
    bait_and_trim(forward_reads='tests/test_fastqs/paired_R1.fastq.gz',
                  plasmid_db='tests/test_fasta/dummy_db.fasta',
                  output_dir='tests/out',
                  low_memory=True)
    assert os.path.isfile('tests/out/plasmid_reads_R1.fastq.gz')
    shutil.rmtree('tests/out')


def test_fasta_write():
    create_individual_fastas(plasmid_db='tests/test_fasta/dummy_db.fasta',
                             potential_plasmid_list=['seq1'],
                             output_dir='tests/fasta/')
    assert os.path.isfile('tests/fasta/seq1') and not os.path.isfile('tests/fasta/seq2')
    shutil.rmtree('tests/fasta')


def test_fasta_kmerization():
    kmerize_individual_fastas(potential_plasmid_list=['dummy_db.fasta'],
                              fasta_dir='tests/test_fasta',
                              output_dir='tests/kmerization')
    assert os.path.isfile('tests/kmerization/dummy_db.fasta.kmc_pre') and os.path.isfile('tests/kmerization/dummy_db.fasta.kmc_suf')
    shutil.rmtree('tests/kmerization')
