from biotools import accessoryfunctions
from biotools import kmc
import multiprocessing
import glob
import os


def compare(plasmid, read_db):
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
    if percent > 0.99:
        return plasmid
    else:
        return 'NA'


def find_plasmids(reads, kmer_db, threads=4):
    plasmids = list()
    kmc.kmc(reads, 'read_db', min_occurrences=2)
    plasmid_files = glob.glob(kmer_db + '/*.kmc_pre')
    read_list = ['read_db'] * len(plasmid_files)
    pool = multiprocessing.Pool(processes=threads)
    results = pool.starmap(compare, zip(plasmid_files, read_list))
    pool.close()
    pool.join()
    for result in results:
        if result != 'NA':
            plasmids.append(result)
    return plasmids
