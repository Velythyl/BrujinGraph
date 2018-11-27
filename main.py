"""
# test de colisions avec notre fonction de hashage
import itertools
alphabet = ['A', 'T', 'C', 'G']
str_tab = []
for i in itertools.product(alphabet, repeat = 6):
    str_tab.append(''.join(i))
print(len(str_tab))
from HashTab import HashTabADN
test_set = str_tab[:2000]
temp = HashTabADN(test_set, 6)
temp.colision_test(test_set)
"""
# Function vraiment pas efficace: je le sais. Elle n'a besoin que de rouler une seule fois pour trouver les primes entre
# 0 et 4^21 pour notre fonction de compression. Après, la liste de primes est mise en fichier.
import gzip
import time
import multiprocessing as mp
from BrujinGraph import DeBrujinGraph

use_mp = True

def read_fasta(path):
    with gzip.open(path, 'rt') as f:
        accession, description, seq = None, None, None
        for line in f:
            if line[0] == '>':
                # yield current record
                if accession is not None:
                    yield seq
                    # break #TEMP TODO
                # start a new record
                accession, description = line[1:].rstrip().split(maxsplit=1)
                seq = ''
            else:
                seq += line.rstrip()


if __name__ == '__main__':
    test = list(read_fasta('GCF_000002985.6_WBcel235_rna.fna.gz'))
    print("FASTA loaded. Building graph...")

    start = time.time()
    test = DeBrujinGraph(test, use_mp=use_mp)
    print("Graph built in", time.time() - start, "seconds.")
