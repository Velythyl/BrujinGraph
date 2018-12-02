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
import gc
import gzip
import time
import multiprocessing as mp
import BrujinGraph
from BrujinGraph import DeBrujinGraph

use_mp = True


def read_fasta(path, only_seq=False):
    with gzip.open(path, 'rt') as f:
        accession, description, seq = None, None, None
        for line in f:
            if line[0] == '>':
                # yield current record
                if accession is not None:
                    if only_seq:
                        yield seq
                    else:
                        yield accession, description, seq
                # start a new record
                accession, description = line[1:].rstrip().split(maxsplit=1)
                seq = ''
            else:
                seq += line.rstrip()

counter = 0
def get_all_kmer(path):
    global counter

    with gzip.open(path, 'rt') as f:
        for line in f:
            seqid, description = line[1:].rstrip().split(maxsplit=1)

            sequence = f.readline().rstrip()
            _ = f.readline()
            quality = f.readline().rstrip()
            for kmer in BrujinGraph.build_kmers(sequence):
                yield kmer
            print("\t", counter, "kmerised.")
            counter += 1

def read_fastq(path, only_seq=False):
    global counter_r

    with gzip.open(path, 'rt') as f:
        for line in f:
            seqid, description = line[1:].rstrip().split(maxsplit=1)

            sequence = f.readline().rstrip()
            _ = f.readline()
            quality = f.readline().rstrip()
            if only_seq:
                yield sequence
            else:
                yield seqid, description, sequence, quality


def _mp_init(semaphore_):
    global semaphore
    semaphore = semaphore_


def _mp_hash_mapper(node):
    semaphore.acquire()
    kmer_list = []

    for kmer in BrujinGraph.build_kmers(node):
        kmer_list.append(BrujinGraph.hash(kmer))

    return kmer_list

#606 avec mp 763 sans


def _mp_hash_all(graph, iter, pool, semaphore):  # hash est op la plus couteuse, on la
    counter = 0

    iterator = pool.imap_unordered(_mp_hash_mapper, iter, chunksize=500)
    pool.close()
    for node_list in iterator:
        graph.add_node_list(node_list)
        print("\t", counter, "added")
        counter += 1
        semaphore.release()

    pool.join()
    gc.collect()


def print_alpha_numera():
    alpha = ['A', 'T', 'C', 'G']
    numera = ['0', '1', '2', '3']

    for i in range(4):
        for j in range(4):
            for k in range(4):
                mot = alpha[i] + alpha[j] + alpha[k]
                num = numera[i] + numera[j] + numera[k]
                print("'" + mot + "': '" + num + "', '" + num + "': '" + mot + "',")

print(len(list(read_fastq('reads.fastq.gz'))))

if __name__ == '__main__':
    start = time.time()

    # print(sys.getsizeof(786945046841))
    # print(sys.getsizeof("ATGCGAGTCTCCACGTCAGTC"))

    graph = None
    if use_mp:
        # Pourquoi est-ce que les methodes de mp sont dans main.py?
        #
        # Python gere mal le multiprocessing. Lorsque les threads de mp requierent des methodes provenant d'une instance
        # il y a une explosion de memoire. Ceci est arrangeable avec des methodes statiques ou des fonctions, mais en
        # instanciant le pool de mp dans le DeBrujinGraph, on gardait quand meme une copie en memoire de l'instance dans
        # chaque process. Donc, on les instancie ici.
        #
        # Cependant, on a besoin d'un sephamore: meme si hash est la methode la plus lente, on se retrouve avec trop de
        # kmer hashes pour la methode add_node_list: il y a explosion de memoire ici aussi. On ne peut pas instancier
        # le pool ici avec son sephamore et utiliser ce dernier dans DeBrujinGraph: python dit qu'il ne reconnait pas le
        # "name sephamore"... Le plus simple est donc de hasher ici et de passer le resultat au graph, meme si c'est un
        # peu moins propre.
        #
        # Logique du mp: le pool hashe les kmer des noms et retourne une liste de hash par nom. Le process main itere
        # sur ces hash et les ajoute au graph avec add_node_list, qui teste tous les nodes et les ajoute au hash_table
        # ssi le hash n'y est pas deja

        cpu_nb = min(mp.cpu_count(), 3)  # On veut un core juste pour add_node_list (donc qui traite le resultat
        # des autres), et max 3 (car teste avec 3 pour sephamore=2000)

        # https://stackoverflow.com/questions/40922526/memory-usage-steadily-growing-for-multiprocessing-pool-imap-unordered
        lock = mp.Semaphore(2000)
        pool = mp.Pool(cpu_nb, initializer=_mp_init, initargs=(lock,))

        graph = DeBrujinGraph()

        _mp_hash_all(graph, read_fastq('reads.fastq.gz', True), pool, lock)

    else:
        print("FASTA loaded. Building graph...")
        graph = DeBrujinGraph(get_all_kmer('reads.fastq.gz'))

    print("Graph built in", time.time() - start, "seconds.")

    graph.save()

    for ele in graph.walk():
        print(ele)