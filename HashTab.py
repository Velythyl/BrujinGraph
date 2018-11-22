import pickle
import random
from bisect import bisect_left
from BrujinGraph import Node

def _adn_to_int(char):
    if char == "A":
        return 0
    if char == "T":
        return 1
    if char == "C":
        return 2
    if char == "G":
        return 3

class HashTabADNSeed:

    def __init__(self):
        table = []
        self.size = 5120
        load = 0.0


    def hash(self, string):
        key = 0
        for char in string:
            random.seed(key + _adn_to_int(char))
            key = random.randint(0, self.size)
        return key

    def colision_test(self, str_tab):
        dict = {}
        colisions = 0
        for str in str_tab:
            key = self.hash(str)
            try:
                temp = dict[key]
                colisions += 1
            except KeyError:
                dict[key] = 0

        print(colisions)
        return colisions

class HashTabADNSeed2:

    def __init__(self):
        table = []
        self.size = 5120
        load = 0.0


    def hash(self, string):
        key = 0
        for i in range(1):
            char4 = string[-4:]
            string = string[:-4]
            char4 = self._add_adn(char4)
            random.seed(key + char4)
            key = random.randint(0, self.size)
        random.seed(key + _adn_to_int(self._add_adn(string)))
        key = random.randint(0, self.size)
        return key

    def _add_adn(self, char4):
        added = 0
        for char in char4:
            added += _adn_to_int(char)
        return added

    def colision_test(self, str_tab):
        dict = {}
        colisions = 0
        for str in str_tab:
            key = self.hash(str)
            try:
                temp = dict[key]
                colisions += 1
            except KeyError:
                dict[key] = 0

        print(colisions)
        return colisions

class AlphabetError(Exception):
    pass

class ADNHashError(Exception):
    pass

class ADNCompressionError(Exception):
    pass

class Bucket:  # Creer qui si colision, sinon juste val direct (reduit besoin mem)
    def __init__(self, val):
        self._list = []
        self._add(val)

    def add(self, val):
        self._list.append(val)

    def __iter__(self):
        for ele in self._list:
            yield ele

    def __len__(self):
        return len(self._list)


class HashTabADN:

    def __init__(self, iter = None, size = 1000, word_length = 21, factor = 2):
        if iter is None:
            self._size = size * factor  # length iterable * 1.25...
            self._used = 0
        else:
            self._size = len(iter)*factor
            self._used = len(iter)

        self._table = [None]*self._size

        self.load = self._used / self._size
        self._safe_bound = 4**word_length       # Nb total de mots dans notre langage

        self._prime = self._get_lowest_prime()
        self._scale = self._get_scale()

        if iter is not None:
            for item in iter:
                self.add(item, True)

    # Pour des questions d'optimisation du resize, on hash les nodes dans le BrujinGraph, a leur creation. De cette
    # facon, toutes les nodes sont deja hashees et on n'a qu'a les compresser pour les entrer dans la table, peu
    # importe si elles sont ajoutees pour la premiere fois a la table ou non
    def add(self, node, is_init = False):
        hashed = node.hash

        com_hash = self._compress(hashed)
        old = self._table[com_hash]

        if isinstance(old, Node):
            bucket = Bucket(old)
            bucket.add(node)
        elif isinstance(old, Bucket):
            old.add(node)
        else:
            self._table[com_hash] = node

        if not is_init:
            self._used += 1
            self._resize()

    def remove(self, name):
        node = self.search(name)

        node_to_refactor = None                 # pour garder optimization pas bucket partout
        com_hash = self._hash_and_compress(name)
        temp_bucket = self._table[com_hash]
        if isinstance(temp_bucket, Bucket):
            temp_bucket._list.remove(node)

            if len(temp_bucket) == 1:
                node_to_refactor = temp_bucket._list[0]
                temp_bucket._list.remove(node_to_refactor)

        self._table[com_hash] = node_to_refactor    # None si pas besoin de refactor donc delete node, et sinon
                                                    # remplace node :)

        self._used -= 1
        self._resize()

    def _resize(self):
        self.load = self._used / self._size
        if self.load >= 0.75:
            self = HashTabADN(self)

    def __len__(self):
        return self._size

    def __iter__(self):
        for ele in self._table:
            if ele is None:
                continue
            elif isinstance(ele, Node):
                yield ele
            elif isinstance(ele, Bucket):
                for sub_ele in ele:
                    yield sub_ele
            else:
                raise Exception

    def search(self, name) -> Node:
        hashed = self.hash(name)                    # prend hash non compresse
        item = self._table[self._compress(hashed)]  # on prend l'item correspondant dans la table
        if item is None:                            # Si None, erreur
            raise KeyError

        if isinstance(item, Node) and item.hash == hashed:                  # Si est une Node, facile
            return item
        elif isinstance(item, Bucket):              # Si est Bucket,
            for ele in item:                        # Cherche dans le bucket pour bonne clef
                if ele.hash == hashed:
                    return ele

        raise KeyError

    """"""
    def _get_lowest_prime(self):    # prend prime approprie selon le safe bound. Peu efficace, mais la recherche ne
        # prends pas longtemps normalement
        for candidat in range(self._size, self._safe_bound):
            is_prime = True
            for facteur in range(2, candidat):
                if candidat % facteur == 0:
                    is_prime = False
                    break
            if is_prime:
                return candidat
        return max(self._size, self._safe_bound)

    def _get_scale(self):
        low = self._prime // 2
        high = low + 1
        while True:
            if low % self._prime != 0:
                return low
            if high % self._prime != 0:
                return high

            low -= 1
            high +=1

            if high >= self._prime:
                raise ADNCompressionError

    """key = 0

        for index in range(len(string)):
            key += (_adn_to_int(string[-index]) * (10 ** index))

        try:
            key = int(str(key), 4)
        except Exception:
            raise ADNHashError
            """
    def hash(self, string):
        key = ""
        for char in string:
            key += self._adn_to_char_int(char)

        try:
            key = int(key, 4)   # base 4! 4 letrres!
        except Exception:
            raise ADNHashError

        return key

    def _adn_to_char_int(self, char):
        if char == "A":
            return '0'
        if char == "T":
            return '1'
        if char == "C":
            return '2'
        if char == "G":
            return '3'
        raise AlphabetError

    def _compress(self, key):
        if self._safe_bound > self._prime:
            key = key*self._scale       # M
            key = key % self._prime     # D
        return key % self._size         # Resize

    def _hash_and_compress(self, string):
        return self._compress(self.hash(string))

    def colision_test(self, str_tab):
        testSet = set()
        colisions = 0
        for str in str_tab:
            key = self._hash_and_compress(str)
            if key in testSet:
                colisions += 1
            else:
                testSet.add(key)

        print("Colisions:",colisions)
        return colisions