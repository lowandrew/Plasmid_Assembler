import glob
import copy
import multiprocessing
import shutil
import argparse
import subprocess
import os
from biotools import kmc
from biotools import accessoryfunctions


class PlasmidScore:
    def __init__(self, plasmid_info):
        self.plasmid = plasmid_info[0]
        self.score = plasmid_info[1]


# This is PlasmidExtractor. It will get plasmid sequences out of raw reads.
def compare(plasmid, read_db, cutoff):
    plasmid = plasmid.replace('.kmc_pre', '')
    db = os.path.split(plasmid)[-1]
    # kmc.kmc(plasmid, db, fm='', tmpdir=db + 'tmp')
    kmc.intersect(plasmid, read_db, db + '_intersect')
    kmc.dump(plasmid, db + '_dump')
    kmc.dump(db + '_intersect', db + '_intersect_dump')
    plasmid_kmers = accessoryfunctions.file_len(db + '_dump')
    read_kmers = accessoryfunctions.file_len(db + '_intersect_dump')
    percent = float(read_kmers)/float(plasmid_kmers)
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


def find_plasmids(reads, kmer_db, cutoff, threads=4):
    plasmids = list()
    kmc.kmc(reads, 'read_db', min_occurrences=2)
    plasmid_files = glob.glob(kmer_db + '/*.kmc_pre')
    read_list = ['read_db'] * len(plasmid_files)
    cutoff_list = [cutoff] * len(plasmid_files)
    pool = multiprocessing.Pool(processes=threads)
    results = pool.starmap(compare, zip(plasmid_files, read_list, cutoff_list))
    pool.close()
    pool.join()
    for result in results:
        if result[0] != 'NA':
            x = PlasmidScore(result)
            plasmids.append(result[0])
            # print(x.plasmid)
            # print(x.score)
    return plasmids


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


def filter_plasmids(plasmids, tmpdir, kmer_db, sequence_db, logfile='logfile.log'):
    """
    When given a list of plasmid files, figures out which ones are very similar to each other and can be ignored,
    picking only one of them.
    :param plasmids: List of paths to plasmid files (these must be in fasta format!).
    :param tmpdir: Temporary directory where some stuff will get placed.
    :param kmer_db: Directory with kmerized plasmid files.
    :param sequence_db: Directory with sequence files in fasta format.
    :param logfile: Logfile to write stdout and stderr of commands to.
    :return: List of plasmids that are relatively unrelated.
    """
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


def check_dependencies():
    """
    Makes sure all the programs that need to be installed are installed.
    """
    dependencies = ['samtools', 'bcftools', 'bedtools', 'bbmap.sh', 'bbduk.sh', 'mash', 'kmc']
    for dep in dependencies:
        is_present = shutil.which(dep)
        if is_present is None:
            raise ModuleNotFoundError('ERROR! Could not find executable for: {}!'.format(dep))


class PlasmidExtractor:
    def __init__(self, reads, args):
        self.forward_reads = reads[0]
        self.reverse_reads = reads[1]
        self.output_base = os.path.join(args.output_dir, os.path.split(self.forward_reads)[-1].split('_')[0])
        self.threads = args.threads
        self.keep_tmpfiles = args.keep_tmpfiles
        self.logfile = self.output_base + '.log'
        self.kmer_db = args.kmer_db
        self.sequence_db = args.sequence_db
        self.cutoff = args.cutoff

    def main(self):
        self.only_reads()
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
        plasmids = find_plasmids(tmpdir + 'interleaved.fastq', self.kmer_db, self.cutoff, threads=self.threads)
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
    args = parser.parse_args()
    paired_reads = find_paired_reads(args.input_directory)
    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir)
    for pair in paired_reads:
        extractor = PlasmidExtractor(pair, args)
        extractor.main()
    # Now need to do some post-analysis on the plasmids found (optionally) to figure out things like how similar they
    # are, what resistance/virulence genes they have, and all that fun stuff.
