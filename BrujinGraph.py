import gc
import math
import os
import pickle
import sys
import time
import sys

sys.setrecursionlimit(10000)


# Opere sur des noms/kmer/str

class DeBrujinGraph:
    def __init__(self, nodes=None, k=21):
        self._k = k
        self.hash_table = HashTabADN()
        # self.hash_table.colision_test(nodes)
        if nodes is not None:
            for kmer in nodes:
                self.add(kmer)

    def __contains__(self, N: str) -> bool:  # détermine si le graphe de Brujin contient le noeud N
        try:
            self.hash_table.search_str(N)
            return True
        except KeyError:
            return False

    def __iter__(self):
        return self.nodes()  # retourne un itérable sur les noeuds du graphe

    def load_factor(self) -> float:  # calcule le facteur de charge de la table de hachage sous-jacente
        return self.hash_table.load

    def add(self, N: str):  # ajoute le noeud N au graphe
        node = hash(N)
        if node not in self.hash_table:
            self.hash_table.add(node)

    def add_node_list(self, node_list):  # ajoute le noeud N au graphe
        for node in node_list:
            if node not in self.hash_table:
                self.hash_table.add(node)

    def _get_succ(self, kmer: str):
        name = kmer[1:]
        for char in ['A', 'T', 'C', 'G']:
            yield name + char

    def _get_pred(self, kmer: str):
        name = kmer[:-1]
        for char in ['A', 'T', 'C', 'G']:
            yield char + name

    def remove(self, N: str):  # enlève le noeud N du graphe
        self.hash_table.remove(self.hash(N))

    def nodes(self):  # retourne un itérable sur les noeuds du graphe
        for node in self.hash_table:
            yield unhash(node)

    def predecessors(self, N: str):  # retourne tous les prédécesseur du noeud N
        return self.cessors(N, True)

    def successors(self, N: str):  # retourne tous les successeurs du noeud N
        return self.cessors(N, False)

    def cessors(self, kmer: str, is_pred: bool):
        cessor_list = []
        for item in self._get_pred(kmer) if is_pred else self._get_succ(kmer):
            try:
                self.hash_table.search_str(item)
                cessor_list.append(item)
            except KeyError:
                pass
        return cessor_list

    def walk(self):
        closed = set()
        while len(closed) != self.hash_table._used:
            start = None
            for kmer in self.hash_table:
                temp = unhash(kmer)
                if (len(self.predecessors(temp)) == 0) and (temp not in closed):
                    start = temp
                    break

            if start == None:
                break   # TODO

            temp_set = set()
            temp_set.add(start)
            for tuple in self._walk(start, temp_set, start):
                for ele in tuple[1]:
                    closed.add(ele)
                yield tuple[0]

    def _walk(self, at, closed, contig):
        succ_list = self.successors(at)
        if len(succ_list) == 0:
            yield (contig, closed)
        else:
            for succ in self.successors(at):
                if succ in closed:
                    yield (contig + succ[-1], closed)
                else:
                    closed.add(succ)
                    # print(((len(contig) - len(succ)) * ' ') + succ)
                    yield from self._walk(succ, closed, contig + succ[-1])

    def save(self, f="DBG.gra"):
        with open(f, "wb") as file:
            pickle.dump(self.__dict__, file)

    @staticmethod
    def loader(f="DBG.gra"):
        with open(f, "rb") as file:
            return pickle.load(file)


def build_kmers(name, k=21):
    for i in range(len(name) - k + 1):
        yield name[i:i + k]


# Opere sur des int/hash/compressions

class ADNCompressionError(object):
    pass

class HashTabADN:

    def __init__(self, size=400000, word_length=21):
        print("\tInitial creation of table:")

        self._size = size  # length iterable * 2...
        self._used = 0

        print("\t\t\tsize: ", self._size)

        self._table = [None] * self._size  # False: rien n'a ete la, True: l'endroit est libre mais etait pris avant

        self.load = self._used / self._size
        self._safe_bound = 4 ** word_length  # Nb total de mots dans notre langage

        self._prime = 92821  # self._get_lowest_prime()
        self._scale = 46410  # self._get_scale()

    def copy(self):
        new = HashTabADN()
        new.__dict__.update(self.__dict__)
        return new

    def _rebuild(self, old, old_used, factor=4):
        print("\t\tResized: rebuilding...")

        self._size = old_used * factor
        print("\t\t\tNew size:", self._size)
        self._used = old_used
        self._table = [None] * self._size

        for item in old:
            self.add(item, True)

    # Pour des questions d'optimisation du resize, on hash les nodes dans le BrujinGraph, a leur creation. De cette
    # facon, toutes les nodes sont deja hashees et on n'a qu'a les compresser pour les entrer dans la table, peu
    # importe si elles sont ajoutees pour la premiere fois a la table ou non
    def add(self, node, is_init=False):

        com_hash = self._compress(node)
        old = self._table[com_hash]

        if old is None:
            self._table[com_hash] = node
        elif isinstance(old, list):
            old.append(node)
        else:
            bucket = [old]
            bucket.append(node)
            self._table[com_hash] = bucket

        if not is_init:
            self._used += 1
            self._resize()

    def remove(self, node):
        place = self.search_node(node)  # lance KeyError si node existe pas

        node_to_refactor = None  # pour garder optimization pas bucket partout
        com_hash = self._compress(node)

        temp_bucket = self._table[com_hash]
        if isinstance(temp_bucket, list):
            temp_bucket.remove(node)

            if len(temp_bucket) == 1:
                node_to_refactor = temp_bucket[0]
            else:
                return

        self._table[com_hash] = node_to_refactor  # None si pas besoin de refactor donc delete node, et sinon
        # remplace node!

        self._used -= 1

    def remove_str(self, str):
        hashed = hash(str)
        if hashed in self:
            self.remove(hashed)

    def _resize(self):
        self.load = self._used / self._size
        if self.load >= 0.75:
            self._rebuild(self, self._used)

    def __len__(self):
        return self._used

    def __contains__(self, node):
        try:
            self.search_node(node)
            return True
        except KeyError:
            return False

    def size(self):
        return self._size

    def __iter__(self):
        for ele in self._table:
            try:
                for ele2 in ele:
                    yield ele2
            except Exception:
                if ele is None:
                    continue
                yield ele

    def search_node(self, node) -> int:
        com_hash = self._compress(node)
        item = self._table[com_hash]  # on prend l'item correspondant dans la table

        if item == node:  # Si est une Node, on teste si c'est la meme
            return com_hash

        try:
            if node in item:
                return node
        except Exception:
            pass

        raise KeyError

    def search_str(self, name) -> int:  # Retourne le hash du str SI il existe dans la table
        hashed = hash(name)  # prend hash non compresse
        return self.search_node(hashed)

    # Prend prime approprie selon le safe bound. Peu efficace, mais rend notre compression
    # bien meilleure.
    #
    # On ne fait pas de cas speciaux pour les petits primes evidents car on commence toujours avec une borne inferieure
    # de 1000 (grosseur de table minimum) de toute facon
    def _get_lowest_prime(self):
        for candidat in range(self.size(), self._safe_bound):
            i = 2
            is_prime = True
            while i * i <= candidat:
                if candidat % i == 0:
                    is_prime = False
                i += 1
            if is_prime:
                print("\t\t\tprime: ", candidat)
                return candidat

    def _get_scale(self):
        low = self._prime // 2
        high = low + 1
        while high <= self._prime:
            if low % self._prime != 0:
                print("\t\t\tscale: ", low)
                return low
            if high % self._prime != 0:
                print("\t\t\tscale: ", high)
                return high

            low -= 1
            high += 1

        raise ADNCompressionError

    def _compress(self, key):
        if self._safe_bound > self._prime:
            key = key * self._scale  # M
            key = key % self._prime  # D
        return key % self._size  # Resize

    def _hash_and_compress(self, string):
        return self._compress(hash(string))

    def colision_test(self, str_tab):
        test_set = set()
        colisions = 0
        for string in str_tab:
            key = self._hash_and_compress(string)
            if key in test_set:
                colisions += 1
            else:
                test_set.add(key)

        print("\tColisions:", colisions)
        return colisions


conv_dict = {'A': '0', 'C': '1', 'T': '2', 'G': '3', '0': 'A', '1': 'C', '2': 'T', '3': 'G'}


def unhash(key):
    # loosely inspire de https://www.codevscolor.com/python-convert-decimal-ternarybase-3/
    def _to_quaternary(num):  # 2
        q, r = divmod(num, 4)
        if q == 0:  # 4
            return str(r)
        else:
            return _to_quaternary(q) + str(r)

    key = _to_quaternary(key)
    string = ""
    for char in key:
        string += conv_dict[char]

    while (len(string) < 21):
        string = 'A' + string

    return string


class ADNHashError(Exception):
    pass


def hash(string):
    try:
        key = ""
        for char in string:
            key += conv_dict[char]

        key = int(key, 4)  # base 4! 4 letrres!

        return key
    except Exception:
        raise ADNHashError