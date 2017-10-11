import bisect
import multiprocessing
import glob
# Attempt to compare two kmer count files, and see how often the kmers from one are in another.


def compare(plasmid_file, read_kmers):
    plasmid_kmers = list()
    f = open(plasmid_file)
    lines = f.readlines()
    f.close()
    # plasmid_kmers = np.empty(len(lines), dtype='U40')
    for line in lines:
        # bisect.insort(plasmid_kmers, line.split()[0])
        plasmid_kmers.append(line.split()[0])
    plasmid_kmers = sorted(plasmid_kmers)
    num_plasmid_kmers = len(plasmid_kmers)
    count = 0
    for kmer in read_kmers:
        idx = bisect.bisect(plasmid_kmers, kmer)
        if plasmid_kmers[idx - 1] == kmer:
            count += 1

    percentage_hit = float(count)/float(num_plasmid_kmers)
    if percentage_hit > 0.8:
        return plasmid_file
    else:
        return 'NA'


def get_read_kmer_list(read_file):
    read_kmers = list()
    f = open(read_file)
    lines = f.readlines()
    f.close()
    for line in lines:
        # bisect.insort(read_kmers, line.split()[0])
        read_kmers.append(line.split()[0])
    read_kmers = sorted(read_kmers)
    return read_kmers


def find_plasmids(read_counts):
    plasmids = list()
    read_kmers = get_read_kmer_list(read_counts)
    plasmid_files = glob.glob('kmerized_plasmids/*')
    read_list = [read_kmers] * len(plasmid_files)
    pool = multiprocessing.Pool(processes=12)
    results = pool.starmap(compare, zip(plasmid_files, read_list))
    pool.close()
    pool.join()
    for result in results:
        if result != 'NA':
            plasmids.append(result)
    return plasmids
    # compare('NZ_LN890526_counts.tab', read_kmers)
