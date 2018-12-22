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
import multiprocessing as mp
import time

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


# Cette fonction initialise le pool de mp pour le hashage: on a besoin d'un semaphore pour se limiter contre main
def _mp_init(semaphore_):
    global semaphore
    semaphore = semaphore_


# Cette fonction hash tous les kmers pouvant etre tires d'un string
def _mp_hash_mapper(string):
    semaphore.acquire()
    kmer_list = []

    for kmer in BrujinGraph.build_kmers(string):
        kmer_list.append(BrujinGraph.hash(kmer))

    return kmer_list


# Cette fonction initialise le pool qui parcours le graph: elle a besoin du graph
def _mp_init_walk(graph_):
    global graph
    graph = graph_


# Pourquoi ne pas simplement retourner temp ssi il n'a pas de predecesseurs? Etrangement, les fonctions de mappage de mp
# Retournent TOUJOURS un element, meme si on n'en retourne pas. De cette facon-ci, si temp_start n'est pas vraiment un
# start, on retourne une liste vide! Comme ca, lorsqu'on construira la liste de starts, on aura simplement rien dans la
# liste.
def _mp_get_starts(ele):
    temp_list = []

    temp_start = BrujinGraph.unhash(ele)
    if len(graph.predecessors(temp_start)) == 0:
        temp_list.append(temp_start)

    return temp_list


# Cette fonction parcours le graph a partir d'un point start
def _mp_walk_mapper(start):
    temp_set = set()
    temp_set.add(start)

    contig_list = []

    for contig in graph._walk(start, temp_set, start):
        contig_list.append(contig)

    return contig_list


# Cette fonction append un fichier avec une "batch" de contigs
def contig_to_file(batch):
    with open("contigs.fa", "a") as file:
        file.write(batch)


use_mp = True

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
        # kmer hashes pour la methode add_hashed_node_list: il y a explosion de memoire ici aussi. On ne peut pas instancier
        # le pool ici avec son sephamore et utiliser ce dernier dans DeBrujinGraph: python dit qu'il ne reconnait pas le
        # "name sephamore"... Le plus simple est donc de hasher ici et de passer le resultat au graph, meme si c'est un
        # peu moins propre.
        #
        # Logique du mp: le pool hashe les kmer des noms et retourne une liste de hash par nom. Le process main itere
        # sur ces hash et les ajoute au graph avec add_hashed_node_list, qui teste tous les nodes et les ajoute au hash_table
        # ssi le hash n'y est pas deja

        cpu_nb = min(mp.cpu_count(), 3)  # On veut un core juste pour add_hashed_node_list (donc qui traite le resultat
        # des autres), et max 3 (car teste avec 3 pour sephamore=2000)

        # https://stackoverflow.com/questions/40922526/memory-usage-steadily-growing-for-multiprocessing-pool-imap-unordered
        semaphore = mp.Semaphore(2000)
        pool = mp.Pool(cpu_nb, initializer=_mp_init, initargs=(semaphore,))

        graph = DeBrujinGraph()

        counter = 0
        iterator = pool.imap_unordered(_mp_hash_mapper, read_fastq('reads.fastq.gz', True), chunksize=500)
        pool.close()
        for node_list in iterator:
            graph.add_hashed_node_list(node_list)
            print("\t", counter, "added")
            counter += 1
            semaphore.release()
            #if counter == 50000: break  # TODO TEMP

        pool.terminate()

    else:
        print("FASTA loaded. Building graph...")
        graph = DeBrujinGraph(get_all_kmer('reads.fastq.gz'))

    print("Graph built in", time.time() - start, "seconds.")

    # graph.save()

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

    if use_mp:
        start = time.time()

        # Ici, c'est le contraire de lorsqu'on utilise le mp pour le hashing: on VEUT que le graph se partage
        # On pourrait donc s'attendre a ce qu'en placant les methodes dans DeBrujinGraph on obitenne le resultat
        # escompte.
        #
        # Oui, le graph se partage de cette facon, mais il doit se re-partager periodiquement, produisant BEAUCOUP
        # d'overhead. En placant le graph comme parametre d'initialisation (dans initargs), on ne fait le partage qu'une
        # seule fois et tout fonctionne. Les methodes de cette mp sont donc aussi dans main pour eviter le overhead

        cpu_nb = min(mp.cpu_count(), 3)

        pool = mp.Pool(cpu_nb, initializer=_mp_init_walk, initargs=(graph,))

        start_list = []
        for temp_list in pool.imap_unordered(_mp_get_starts, graph.hash_table, chunksize=500):
            for temp_start in temp_list:
                start_list.append(temp_start)

        iterator = pool.imap_unordered(_mp_walk_mapper, start_list, chunksize=500)
        pool.close()

        for contig_list in iterator:
            for contig in contig_list:
                # print(contig)
                batch += ">" + build_id() + " DESCRIPTION\n" + contig + "\n"
                counter += 1
                if counter >= batch_size:
                    contig_to_file(batch)
                    batch = ""
                    counter = 0

        pool.terminate()

    else:
        start = time.time()

        for ele in graph.walk():
            # print(ele)
            batch += ">" + build_id() + " DESCRIPTION\n" + ele + "\n"
            counter += 1
            if counter >= batch_size:
                contig_to_file(batch)
                batch = ""
                counter = 0

    if counter < batch_size:
        with open("contigs.fa", "a") as file:
            file.write(batch)

    print("Graph walked in", time.time() - start, "seconds.")
