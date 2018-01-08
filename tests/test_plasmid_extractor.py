import sys
import os
import shutil

parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(0, parentdir)

from New_Extractor import *


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

