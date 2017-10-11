import multiprocessing
import argparse
import sys
import subprocess
import os
import kmer_overlap
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from io import StringIO
from Bio import SeqIO
# This is PlasmidExtractor. It will get plasmid sequences out of raw reads/draft assemblies.

# Things PlasmidExtractor needs to be able to do:
# 1) Figure out which plasmids are present in raw reads/assemblies.
# 2) Map raw reads to closest reference plasmid and generate a consensus sequence.
# 3) Identify contigs that are likely of plasmid origin in a draft assembly.


def mash_fasta(assembly, plasmid_sketch, outdir='tmp/', threads=8):
    """
    :param assembly: Fasta assembly file.
    :param plasmid_sketch: .msh file created by mash sketch
    :param outdir: output directory where tsv file called mash_output containing distances for each contig in assembly
    to each plasmid in the plasmid_sketch file will be stored.
    :param threads: Number of threads to run mash on.
    """
    print('Mashing as initial screen.')
    # Add parallelism eventually.
    cmd = 'mash dist -d 0.1 -p {} -i {} {} > {}/mash_output'.format(assembly, str(threads), plasmid_sketch, outdir)
    print(cmd)
    subprocess.call(cmd, shell=True)


def parse_mash_output(mash_output):
    """
    :param mash_output: Tab-delimited file created by mash.
    :return: dictionary containing best plasmid match to each contig.
    """
    print('Parsing mash output.')
    f = open(mash_output)
    mash_data = f.readlines()
    f.close()
    best_hits = dict()
    distances = dict()
    for match in mash_data:
        x = match.split()
        contig = x[0]
        plasmid = x[1]
        distance = float(x[2])
        if contig not in best_hits:
            best_hits[contig] = plasmid
            distances[contig] = distance
        else:
            if distances[contig] > distance:
                best_hits[contig] = plasmid
                distances[contig] = distance
    return best_hits


def make_blast_database(fasta_file):
    """
    Checks that a fasta file has the files needed for it to be a blast db, and creates the db if it doesn't.
    :param fasta_file: Full path to fasta file.
    :return:
    """
    database_extensions = ['.nhr', '.nin', '.nsq']
    db_exists = True
    for extension in database_extensions:
        if not os.path.isfile(fasta_file + extension):
            db_exists = False
    if not db_exists:
        cmd = 'makeblastdb -dbtype nucl -in {}'.format(fasta_file)
        subprocess.call(cmd, shell=True)


def retrieve_contig_sequence(contig_header, fasta_file):
    """
    :param contig_header: Title for a contig you'd like the sequence of.
    :param fasta_file: Path to fasta file containing the contig.
    :return: NA if contig wasn't found, the sequence of the contig if it was found.
    """
    sequence = 'NA'
    for contig in SeqIO.parse(fasta_file, 'fasta'):
        if contig.id == contig_header:
            sequence = str(contig.seq)
    return sequence


# Name this better.
def blast_things(best_hits, assembly):
    """
    :param best_hits: Dictionary created by parse_mash_output
    :param assembly: Fasta file that has your draft assembly of interest.
    :return: Dictionary containing high-scoring plasmid locations in format:
    {contig_name:[[start1, end1], [start2, end2]]}
    """
    plasmid_locations = dict()
    for contig in best_hits:
        # Check that blast db exists for the plasmid, and create it if it doesn't.
        make_blast_database(best_hits[contig])
        sequence = retrieve_contig_sequence(contig, assembly)
        blastn = NcbiblastnCommandline(db=best_hits[contig], outfmt=5)
        stdout, stderr = blastn(stdin=sequence)
        # print(contig, best_hits[contig])
        locations = list()
        for record in NCBIXML.parse(StringIO(stdout)):
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    # print(hsp.query_start, hsp.query_start + hsp.align_length)
                    locations.append([hsp.query_start, hsp.query_start + hsp.align_length])
                    if hsp.expect > 1e-180:
                        break
        plasmid_locations[contig] = locations
    return plasmid_locations


def quality_trim(forward_in, reverse_in, forward_out, reverse_out):
    """
    Runs bbduk to quality trim a set of reads.
    :param forward_in: Forward reads.
    :param reverse_in: Reverse reads.
    :param forward_out:
    :param reverse_out:
    :return:
    """
    print('Quality Trimming...')
    cmd = 'which bbduk.sh'
    bbduk_dir = subprocess.check_output(cmd.split()).decode('utf-8')
    bbduk_dir = bbduk_dir.split('/')[:-1]
    bbduk_dir = '/'.join(bbduk_dir)
    cmd = 'bbduk.sh in1={} in2={} out1={} out2={} qtrim=w trimq=20 k=25 minlength=50 forcetrimleft=15' \
          ' ref={}/resources/adapters.fa overwrite hdist=1 tpe tbo'.format(forward_in, reverse_in,
                                                                           forward_out, reverse_out, bbduk_dir)
    subprocess.call(cmd, shell=True)


def get_plasmid_reads(forward_in, reverse_in, forward_out, reverse_out, plasmid_database):
    print('Extracting plasmid reads...')
    cmd = 'bbduk.sh in1={} in2={} outm={} outm2={} ref={} overwrite'.format(forward_in, reverse_in, forward_out,
                                                                            reverse_out, plasmid_database)
    subprocess.call(cmd, shell=True)


def interleave_reads(forward_in, reverse_in, interleaved):
    print('Interleaving!')
    cmd = 'reformat.sh in1={} in2={} out={} overwrite'.format(forward_in, reverse_in, interleaved)
    subprocess.call(cmd, shell=True)


def mash_reads(interleaved_reads, plasmid_sketch):
    cmd = 'mash dist -d 0.15 -p 12 -r {} {} > tmp/mash_output'.format(interleaved_reads, plasmid_sketch)
    subprocess.call(cmd, shell=True)


def generate_consensus(reads, reference_fasta, output_fasta):
    # Step 1: Index fasta file
    cmd = 'samtools faidx {}'.format(reference_fasta)
    subprocess.call(cmd, shell=True)
    # Step 2: Run bbmap to generate sam/bamfile.
    cmd = 'bbmap.sh ref={} in={} out={} nodisk overwrite'.format(reference_fasta, reads, 'out.bam')
    subprocess.call(cmd, shell=True)
    # Step 3: Sort the bam file.
    cmd = 'samtools sort out.bam -o out_sorted.bam'
    subprocess.call(cmd, shell=True)
    # Step 3: Fancy bcftools piping to generate vcf file.
    cmd = 'bcftools mpileup -Ou -f {} out_sorted.bam | bcftools call -mv -Oz -o calls.vcf.gz'.format(reference_fasta)
    subprocess.call(cmd, shell=True)
    # Step 4: Index vcf file.
    cmd = 'tabix calls.vcf.gz'
    subprocess.call(cmd, shell=True)
    # Step 5: Generate consensus fasta from vcf file.
    cmd = 'cat {} | bcftools consensus calls.vcf.gz > {}'.format(reference_fasta, output_fasta)
    subprocess.call(cmd, shell=True)


def kmerize_reads(read_file, output_file):
    cmd = 'jellyfish count -m 31 -t 12 -s 100M --bf-size 100M -C ' + read_file
    subprocess.call(cmd, shell=True)
    cmd = 'jellyfish dump -c mer_counts.jf > ' + output_file
    subprocess.call(cmd, shell=True)


def only_assembly(assembly, output_base):
    tmpdir = output_base + 'tmp/'
    if not os.path.isdir(output_base):
        os.makedirs(output_base)
    if not os.path.isdir(tmpdir):
        os.makedirs(tmpdir)
    print('Analyzing draft assembly.')
    # Step 1: Mash each contig against the plasmid database sketch to see what's most closely related and screen out
    # things we don't want to blast against.
    mash_fasta(assembly, 'reduced_sketch.msh', outdir=tmpdir)
    hits = parse_mash_output(tmpdir + '/mash_output')
    # Step 2: Once MASH has done a preliminary screen, get to BLASTing to more precisely locate plasmid sequences.
    locations = blast_things(hits, assembly)
    # Step 3: Create a nice report that tells you everything you need to know.
    f = open(output_base + '/plasmid_locations.csv', 'w')
    f.write('Contig,Plasmid,Start,End\n')
    for contig in hits:
        for location in locations[contig]:
            outstr = '{},{},{},{}\n'.format(contig, hits[contig], str(location[0]), str(location[1]))
            f.write(outstr)
    f.close()


def only_reads(reads, output_base):
    tmpdir = output_base + 'tmp/'
    if not os.path.isdir(output_base):
        os.makedirs(output_base)
    print('Analyzing raw reads!')
    # Step 1: Quality trim reads.
    quality_trim(reads[0], reads[1], tmpdir + 'trimmed_R1.fastq.gz', tmpdir + 'trimmed_R2.fastq.gz')
    # Step 1.5: Bait out reads that match to plasmid stuff.
    get_plasmid_reads(tmpdir + 'trimmed_R1.fastq.gz', tmpdir + 'trimmed_R2.fastq.gz',
                      tmpdir + 'plasmid_R1.fastq.gz', tmpdir + 'plasmid_R2.fastq.gz', 'nucleotideseq.fasta')
    interleave_reads(tmpdir + 'plasmid_R1.fastq.gz', tmpdir + 'plasmid_R2.fastq.gz', tmpdir + 'interleaved.fastq')
    # Step 2: MASH reads against plasmid sketch to figure out which plasmids we should be reference mapping.
    # mash_reads('interleaved.fastq.gz', 'plasmid_sketch.msh')
    kmerize_reads(tmpdir + 'interleaved.fastq', tmpdir + 'reads.tab')
    plasmids = kmer_overlap.find_plasmids(tmpdir + 'reads.tab')
    # plasmids = parse_read_mash_output('tmp/mash_output')
    # Step 3: For each very likely plasmid, do reference mapping and generate a consensus sequence.
    for plasmid in plasmids:
        generate_consensus(tmpdir + 'interleaved.fastq', plasmid.replace('kmerized_plasmids', 'reduced_db'),
                           output_base + '/' + os.path.split(plasmid)[-1] + '_consensus.fasta')


class PlasmidExtractor(object):
    def __init__(self, args):
        self.forward_reads = args.reads[0]
        self.reverse_reads = args.reads[1]
        self.output_base = args.output_dir
        self.threads = args.threads
        self.keep_tmpfiles = args.keep_tmpfiles


if __name__ == '__main__':
    num_cpus = multiprocessing.cpu_count()
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--assembly', default='NA', type=str, help='Full path to draft assembly.')
    parser.add_argument('-r', '--reads', default='NA', nargs='+', type=str, help='Raw read files. If paired, enter '
                                                                                 'both, separated by a space.')
    parser.add_argument('output_dir', type=str, help='Where results will be stored.')
    parser.add_argument('-t', '--threads', default=num_cpus, type=int, help='Number of CPUs to run analysis on.'
                                                                            ' Defaults to number of cores on your'
                                                                            ' machine.')
    parser.add_argument('-k', '--keep-tmpfiles', default=False, action='store_true', help='When specified, will keep'
                                                                                          ' the created tmp directory'
                                                                                          ' instead of deleting it.')
    args = parser.parse_args()
    if args.assembly == 'NA' and args.reads == 'NA':
        print('No assembly files or raw reads specified. Exiting...')
        sys.exit()
    elif args.assembly == 'NA' and args.reads != 'NA':
        only_reads(args.reads, args.output_dir)
    elif args.assembly != 'NA' and args.reads == 'NA':
        only_assembly(args.assembly, args.output_dir)
    else:
        print('Input is both raw reads and assembly. I have not fully figured out what it is that I want to do with '
              'this yet.')
