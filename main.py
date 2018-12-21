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

def get_all_kmer(path):
    counter = 0

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


def _mp_init_hash(semaphore_):
    global semaphore
    semaphore = semaphore_

def _mp_init_walk(graph_):
    global graph
    graph = graph_


def _mp_hash_mapper(node):
    semaphore.acquire()
    kmer_list = []

    for kmer in BrujinGraph.build_kmers(node):
        kmer_list.append(BrujinGraph.Node(BrujinGraph.hash(kmer), kmer))

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
        if counter == 50000: return  # TODO TEMP

    pool.join()
    gc.collect()

def _mp_walk_mapper(start):
    temp_set = set()
    temp_set.add(start)

    contig_list = []

    for contig in graph._walk(start, temp_set, start):
        contig_list.append(contig)

    return contig_list

def contig_to_file(batch):
    with open("contigs.fa", "a") as file:
        file.write(batch)

use_mp = True

# Pour faire le test de la question 2.b, mettre True au lieu de False
if False:
    start = time.time()
    d = {}
    #for ele in get_all_kmer('reads.fastq.gz'):
    #    d[ele] = ele
    print("Standard hash table bulds in:", time.time() - start)

    start = time.time()
    DeBrujinGraph(get_all_kmer('reads.fastq.gz'))
    print("Custom hash table builds in:", time.time() - start)



elif __name__ == '__main__':
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
        pool = mp.Pool(cpu_nb, initializer=_mp_init_hash, initargs=(lock,))

        graph = DeBrujinGraph()

        _mp_hash_all(graph, read_fastq('reads.fastq.gz', True), pool, lock)
        pool.terminate()

    else:
        print("FASTA loaded. Building graph...")
        graph = DeBrujinGraph(get_all_kmer('reads.fastq.gz'))

    print("Graph built in", time.time() - start, "seconds.")

    #graph.save()

    # build_id genere des id uniques
    id_counter = 0
    def build_id():
        global id_counter
        temp = str(id_counter)
        while len(temp) < 10:
            temp = "0" + temp

        id_counter += 1
        return temp


    with open("contigs.fa", "w") as file:
        file.write("")  # wipes previous file

    counter = 0
    batch = ""
    batch_size = 5000



    if True:
        start = time.time()

        cpu_nb = min(mp.cpu_count(), 3)

        pool = mp.Pool(cpu_nb, initializer=_mp_init_walk, initargs=(graph,))

        start_list = graph.get_all_starts()

        iterator = pool.imap_unordered(_mp_walk_mapper, start_list, chunksize=1)
        pool.close()

        for contig_list in iterator:
            for contig in contig_list:
                #print(counter, contig)
                batch += ">" + build_id() + " DESCRIPTION\n" + contig + "\n"
                counter += 1
                if counter >= batch_size:
                    contig_to_file(batch)
                    batch = ""
                    counter = 0

        if counter < batch_size:
            with open("contigs.fa", "a") as file:
                file.write(batch)

        pool.join()
        pool.terminate()

        print(time.time() - start, "seconds using mp")

    if True:
        start = time.time()

        for ele in graph.walk():
            #print(counter, ele)
            batch += ">" + build_id() + " DESCRIPTION\n" + ele + "\n"
            counter += 1
            if counter >= batch_size:
                contig_to_file(batch)
                batch = ""
                counter = 0

        if counter < batch_size:
            with open("contigs.fa", "a") as file:
                file.write(batch)

        print(time.time() - start, "seconds using sp")


