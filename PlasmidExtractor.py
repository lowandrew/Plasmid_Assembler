import copy
import multiprocessing
import shutil
import argparse
import sys
import subprocess
import os
import kmer_overlap
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from io import StringIO
from Bio import SeqIO
import kmc_overlap
# This is PlasmidExtractor. It will get plasmid sequences out of raw reads/draft assemblies.

# Things PlasmidExtractor needs to be able to do:
# 1) Figure out which plasmids are present in raw reads/assemblies.
# 2) Map raw reads to closest reference plasmid and generate a consensus sequence.
# 3) Identify contigs that are likely of plasmid origin in a draft assembly.


def mash_fasta(assembly, plasmid_sketch, outdir='tmp/', threads=8, logfile='logfile.log'):
    """
    :param assembly: Fasta assembly file.
    :param plasmid_sketch: .msh file created by mash sketch
    :param outdir: output directory where tsv file called mash_output containing distances for each contig in assembly
    to each plasmid in the plasmid_sketch file will be stored.
    :param threads: Number of threads to run mash on.
    :param logfile: Log to write stdout and stderr to.
    """
    print('Mashing as initial screen.')
    # Add parallelism eventually.
    cmd = 'mash dist -d 0.1 -p {} -i {} {} > {}/mash_output'.format(str(threads), assembly, plasmid_sketch, outdir)
    with open(logfile, 'a+') as log:
        subprocess.call(cmd, shell=True, stdout=log, stderr=log)


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


def make_blast_database(fasta_file, logfile='logfile.log'):
    """
    Checks that a fasta file has the files needed for it to be a blast db, and creates the db if it doesn't.
    :param fasta_file: Full path to fasta file.
    :param logfile: Log to write stdout and stderr to.
    :return:
    """
    database_extensions = ['.nhr', '.nin', '.nsq']
    db_exists = True
    for extension in database_extensions:
        if not os.path.isfile(fasta_file + extension):
            db_exists = False
    if not db_exists:
        cmd = 'makeblastdb -dbtype nucl -in {}'.format(fasta_file)
        with open(logfile, 'a+') as log:
            subprocess.call(cmd, shell=True, stdout=log, stderr=log)


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


def quality_trim(forward_in, reverse_in, forward_out, reverse_out, logfile='logfile.log', threads=4):
    """
    Runs bbduk to quality trim a set of reads - also removes adapters.
    :param forward_in: Forward reads.
    :param reverse_in: Reverse reads.
    :param forward_out: Outputted forward trimmed reads.
    :param reverse_out: Outputted reverse trimmed reads.
    :param threads: Number of threads to run bbduk with.
    :param logfile: Log to write stdout and stderr to.
    """
    print('Quality Trimming...')
    cmd = 'which bbduk.sh'
    bbduk_dir = subprocess.check_output(cmd.split()).decode('utf-8')
    bbduk_dir = bbduk_dir.split('/')[:-1]
    bbduk_dir = '/'.join(bbduk_dir)
    cmd = 'bbduk.sh in1={} in2={} out1={} out2={} qtrim=w trimq=20 k=25 minlength=50 forcetrimleft=15' \
          ' ref={}/resources/adapters.fa overwrite hdist=1 tpe tbo threads={}'.format(forward_in, reverse_in,
                                                                                      forward_out, reverse_out,
                                                                                      bbduk_dir, str(threads))
    with open(logfile, 'a+') as log:
        subprocess.call(cmd, shell=True, stdout=log, stderr=log)


def get_plasmid_reads(forward_in, reverse_in, forward_out, reverse_out, plasmid_database,
                      logfile='logfile.log', threads=4):
    """
    Use bbduk to extract any reads that match to the plasmid database.
    :param forward_in: Forward input reads.
    :param reverse_in: Reverse input reads.
    :param forward_out: Forward reads that match to the database.
    :param reverse_out: Reverse reads that match to the database.
    :param plasmid_database: Fasta-formatted file with sequences you want things to match to.
    :param threads: Number of threads to run bbduk with.
    :param logfile: Log to write stdout and stderr to.
    """
    print('Extracting plasmid reads...')
    cmd = 'bbduk.sh in1={} in2={} outm={} outm2={} ref={} threads={} overwrite'.format(forward_in, reverse_in,
                                                                                       forward_out, reverse_out,
                                                                                       plasmid_database, str(threads))
    with open(logfile, 'a+') as log:
        subprocess.call(cmd, shell=True, stdout=log, stderr=log)


def interleave_reads(forward_in, reverse_in, interleaved, logfile='logfile.log'):
    """
    Uses bbtools to interleave forward and reverse reads.
    :param forward_in: Forward input reads.
    :param reverse_in: Reverse input reads.
    :param interleaved: Interleaved read file.
    :param logfile: Log to write stdout and stderr to.
    """
    print('Interleaving!')
    cmd = 'reformat.sh in1={} in2={} out={} overwrite'.format(forward_in, reverse_in, interleaved)
    with open(logfile, 'a+') as log:
        subprocess.call(cmd, shell=True, stdout=log, stderr=log)


def generate_consensus(reads, reference_fasta, output_fasta, logfile='logfile.log',
                       cleanup=True, output_base='out', threads=4):
    """
    Generates a consensus fasta give a set of interleaved reads and a reference sequence.
    :param reads: Interleaved reads.
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
        cmd = 'bbmap.sh ref={} in={} out={} nodisk overwrite threads={}'.format(reference_fasta, reads,
                                                                                output_base + '.bam', str(threads))
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


def remove_plasmid_from_reads(plasmid_sequences, forward_in, reverse_in, forward_out, reverse_out,
                              tmpfasta, logfile='logfile.log', threads=4):
    """
    Removes reads that have kmers matching to plasmid sequence.
    :param plasmid_sequences: List of plasmid fasta files. Will be concatenated into 1 fasta.
    :param forward_in: Forward input reads.
    :param reverse_in: Reverse input reads.
    :param forward_out: Forward output reads.
    :param reverse_out: Reverse output reads.
    :param tmpfasta: Where your temporary concatenated fasta will be stored.
    :param logfile: Log to write stdout and stderr to.
    :param threads: Number of threads to run analysis with.
    :return:
    """
    print('Generating plasmid-less reads...')
    with open(tmpfasta, 'w') as outfile:
        for sequence in plasmid_sequences:
            with open(sequence) as infile:
                outfile.write(infile.read())
    cmd = 'bbduk.sh ref={} in={} in2={} out={} out2={} threads={} overwrite'.format(tmpfasta, forward_in, reverse_in,
                                                                                    forward_out, reverse_out,
                                                                                    str(threads))
    with open(logfile, 'a+') as log:
        subprocess.call(cmd, shell=True, stdout=log, stderr=log)


def kmerize_reads(read_file, output_file, kmer_size=31, logfile='logfile.log', threads=4):
    """
    Generates kmers using jellyfish on a read file.
    :param read_file: Read file. Can be single-end or interleaved.
    :param output_file: Name of your output file (will be tab-delimited).
    :param kmer_size: Desired mer size. Default 31.
    :param logfile: Logfile to send stdout and stderr to.
    """
    with open(logfile, 'a+') as log:
        cmd = 'jellyfish count -m {} -t {} -s 100M --bf-size 100M -C {}'.format(str(kmer_size), str(threads), read_file)
        subprocess.call(cmd, shell=True, stdout=log, stderr=log)
        cmd = 'jellyfish dump -c mer_counts.jf > ' + output_file
        subprocess.call(cmd, shell=True, stdout=log, stderr=log)


def filter_plasmids(plasmids, tmpdir, kmer_db, sequence_db, logfile='logfile.log'):
    # Given a list of plasmids, try to figure out which ones are overly similar and should be disregarded.
    print('Filtering plasmids.')
    # Step 1: Sketch all of the plasmids.
    cmd = 'mash sketch -o {}/sketch.msh'.format(tmpdir)
    for plasmid in plasmids:
        cmd += ' ' + plasmid.replace(kmer_db, sequence_db)
    with open(logfile, 'a+') as log:
        subprocess.call(cmd, shell=True, stdout=log, stderr=log)
    # Step 2: For each plasmid, find out if any of the other plasmids are too close.
    keeper_plasmids = copy.deepcopy(plasmids)
    for plasmid in plasmids:
        if plasmid in keeper_plasmids:
            cmd = 'mash dist -d 0.02 {} {}/sketch.msh > {}/distances.tab'.format(plasmid.replace(kmer_db, sequence_db),
                                                                                 tmpdir, tmpdir)
            with open(logfile, 'a+') as log:
                subprocess.call(cmd, shell=True, stdout=log, stderr=log)
            f = open(tmpdir + '/distances.tab')
            lines = f.readlines()
            f.close()
            for line in lines:
                info = line.split()
                if info[0] != info[1]:
                    keeper_plasmids.remove(info[1].replace(sequence_db, kmer_db))
    return keeper_plasmids


class PlasmidExtractor(object):
    def __init__(self, args):
        if args.reads != 'NA':
            self.forward_reads = args.reads[0]
            self.reverse_reads = args.reads[1]
        else:
            self.forward_reads = 'NA'
            self.reverse_reads = 'NA'
        self.output_base = args.output_dir
        self.threads = args.threads
        self.keep_tmpfiles = args.keep_tmpfiles
        self.assembly = args.assembly
        self.logfile = self.output_base + '.log'
        self.kmer_db = args.kmer_db
        self.sequence_db = args.sequence_db

    def main(self):
        if self.assembly == 'NA' and self.forward_reads == 'NA':
            print('No assembly files or raw reads specified. Exiting...')
            sys.exit()
        elif self.assembly == 'NA' and self.forward_reads != 'NA':
            self.only_reads()
        elif self.assembly != 'NA' and self.forward_reads == 'NA':
            self.only_assembly()
        else:
            print('Input is both raw reads and assembly. I have not fully figured out what it is that I want to do with'
                  ' this yet.')

        try:
            if not self.keep_tmpfiles:
                shutil.rmtree(self.output_base + 'tmp')
        except FileNotFoundError:
            pass

    def only_reads(self):
        tmpdir = self.output_base + 'tmp/'
        if not os.path.isdir(self.output_base):
            os.makedirs(self.output_base)
        print('Analyzing raw reads!')
        # Step 1: Quality trim reads.
        quality_trim(self.forward_reads, self.reverse_reads, tmpdir + 'trimmed_R1.fastq.gz',
                     tmpdir + 'trimmed_R2.fastq.gz', logfile=self.logfile, threads=self.threads)
        # Step 1.5: Bait out reads that match to plasmid stuff.
        get_plasmid_reads(tmpdir + 'trimmed_R1.fastq.gz', tmpdir + 'trimmed_R2.fastq.gz',
                          tmpdir + 'plasmid_R1.fastq.gz', tmpdir + 'plasmid_R2.fastq.gz', 'nucleotideseq.fa',
                          logfile=self.logfile, threads=self.threads)
        interleave_reads(tmpdir + 'plasmid_R1.fastq.gz', tmpdir + 'plasmid_R2.fastq.gz', tmpdir + 'interleaved.fastq',
                         logfile=self.logfile)
        # Step 2: MASH reads against plasmid sketch to figure out which plasmids we should be reference mapping.
        # mash_reads('interleaved.fastq.gz', 'plasmid_sketch.msh')
        # kmerize_reads(tmpdir + 'interleaved.fastq', tmpdir + 'reads.tab', logfile=self.logfile, threads=self.threads)
        print('Searching for plasmids in raw reads...')
        # plasmids = kmer_overlap.find_plasmids(tmpdir + 'reads.tab', self.kmer_db, threads=self.threads)
        plasmids = kmc_overlap.find_plasmids(tmpdir + 'interleaved.fastq', self.kmer_db, threads=self.threads)
        # plasmids = parse_read_mash_output('tmp/mash_output')
        # Step 3: For each very likely plasmid, do reference mapping and generate a consensus sequence.
        finished_plasmids = list()
        print('Found {} potential plasmids!'.format(len(plasmids)))
        good_plasmids = filter_plasmids(plasmids, tmpdir, self.kmer_db, self.sequence_db, logfile=self.logfile)
        print('Plasmids filtered! Generating consensus for {} plasmids.'.format(str(len(good_plasmids))))
        for plasmid in good_plasmids:
            print('Generating consensus for {}-like sequence...'.format(plasmid))
            generate_consensus(tmpdir + 'interleaved.fastq', plasmid.replace(self.kmer_db, self.sequence_db),
                               self.output_base + '/' + os.path.split(plasmid)[-1] + '_consensus.fasta',
                               logfile=self.logfile, threads=self.threads)
            finished_plasmids.append(self.output_base + '/' + os.path.split(plasmid)[-1] + '_consensus.fasta')
        # Step 4: Generate raw reads that don't have plasmid reads (ideally).
        out_forward = self.forward_reads.replace('.fastq', '_noplasmid.fastq')
        out_forward = os.path.split(out_forward)[-1]
        out_reverse = self.reverse_reads.replace('.fastq', '_noplasmid.fastq')
        out_reverse = os.path.split(out_reverse)[-1]
        if len(good_plasmids) > 0:
            remove_plasmid_from_reads(finished_plasmids, self.forward_reads, self.reverse_reads,
                                      self.output_base + '/' + out_forward,
                                      self.output_base + '/' + out_reverse,
                                      tmpdir + 'concatenated_plasmid.fasta', threads=self.threads)
        print('Plasmid extraction complete!')

    def only_assembly(self):
        tmpdir = self.output_base + 'tmp/'
        if not os.path.isdir(self.output_base):
            os.makedirs(self.output_base)
        if not os.path.isdir(tmpdir):
            os.makedirs(tmpdir)
        print('Analyzing draft assembly.')
        # Step 1: Mash each contig against the plasmid database sketch to see what's most closely related and screen out
        # things we don't want to blast against.
        mash_fasta(self.assembly, 'reduced_sketch.msh', outdir=tmpdir)
        hits = parse_mash_output(tmpdir + '/mash_output')
        # Step 2: Once MASH has done a preliminary screen, get to BLASTing to more precisely locate plasmid sequences.
        locations = blast_things(hits, self.assembly)
        # Step 3: Create a nice report that tells you everything you need to know.
        f = open(self.output_base + '/plasmid_locations.csv', 'w')
        f.write('Contig,Plasmid,Start,End\n')
        for contig in hits:
            for location in locations[contig]:
                outstr = '{},{},{},{}\n'.format(contig, hits[contig], str(location[0]), str(location[1]))
                f.write(outstr)
        f.close()
        print('Plasmid locating complete!')


if __name__ == '__main__':
    num_cpus = multiprocessing.cpu_count()
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--assembly', default='NA', type=str, help='Full path to draft assembly.')
    parser.add_argument('-r', '--reads', default='NA', nargs='+', type=str, help='Raw read files. If paired, enter '
                                                                                 'both, separated by a space.')
    parser.add_argument('output_dir', type=str, help='Where results will be stored.')
    parser.add_argument('kmer_db', type=str, help='Folder containing kmerized plasmids.')
    parser.add_argument('sequence_db', type=str, help='Folder containing plasmid sequences.')
    parser.add_argument('-t', '--threads', default=num_cpus, type=int, help='Number of CPUs to run analysis on.'
                                                                            ' Defaults to number of cores on your'
                                                                            ' machine.')
    parser.add_argument('-k', '--keep-tmpfiles', default=False, action='store_true', help='When specified, will keep'
                                                                                          ' the created tmp directory'
                                                                                          ' instead of deleting it.')
    args = parser.parse_args()
    extractor = PlasmidExtractor(args)
    extractor.main()

