class Node:

    def __init__(self, string, id):
        self.kmer = string
        self.id = id


class DeBrujinGraph:
    id_counter = 0

    @staticmethod
    def _build_id():
        temp = str(DeBrujinGraph.id_counter)
        while len(temp) < 10:
            temp = "" + temp

        DeBrujinGraph.id_counter += 1
        return temp

    def __init__(self, nodes, k=21):
        self.hash_table = HashTabADN(size=int(len(nodes)*1.35))
        self._alphabet = ['A', 'C', 'T', 'G']
        self._k = k
        #self.hash_table.colision_test(nodes)
        i = 0
        l = len(nodes)
        for name in nodes:
            for kmer in self._build_kmers(name):
                self.add(kmer)
            #print("\t", i, " of ", l, " added.")
            i +=1

    def _build_kmers(self, name):
        for i in range(len(name) - self._k + 1):
            yield name[i:i + self._k]

    def _hash(self, string):  # Sucre syntaxique
        return self.hash_table.hash(string)

    def __contains__(self, name: str) -> bool:  # détermine si le graphe de Brujin contient le noeud N
        try:
            self.hash_table.search(name)
            return True
        except KeyError:
            return False

    def __iter__(self):
        return self.nodes()  # retourne un itérable sur les noeuds du graphe TODO NOEUD OU STR???

    def load_factor(self) -> float:  # calcule le facteur de charge de la table de hachage sous-jacente
        return self.hash_table.load

    def add(self, kmer: str):  # ajoute le noeud N au graphe   TODO NOEUD OU STR

        node = Node(kmer, self.hash_table.hash(kmer))
        self.hash_table.add(node)
        """
        for candidate in self._get_pred(kmer):
            try:  # Duck typing
                pred = self.hash_table.search(candidate[0])
                Node.link(pred, node)
            except KeyError:
                pass

        for candidate in self._get_succ(kmer):
            try:  # Duck typing
                succ = self.hash_table.search(candidate[0])
                Node.link(node, succ)
            except KeyError:
                pass
        """

    def _get_succ(self, kmer: str):
        name = kmer[1:]
        for char in self._alphabet:
            yield [name + char, char]

    def _get_pred(self, kmer: str):
        name = kmer[:-1]
        for char in self._alphabet:
            yield [char + name, char]

    def remove(self, kmer: str):  # enlève le noeud N du graphe
        self.hash_table.remove(kmer)

    def nodes(self):  # retourne un itérable sur les noeuds du graphe TODO NOEUD OU STR???
        for node in self.hash_table:
            yield node.name

    def predecessors(self, kmer: str):  # retourne tous les prédécesseur du noeud N
        return self.cessors(kmer, True)

    def successors(self, kmer: str):  # retourne tous les successeurs du noeud N
        return self.cessors(kmer, False)

    def cessors(self, kmer: str, is_pred: bool):
        cessor_list = []
        for pred in self._get_pred(kmer) if is_pred else self._get_succ(kmer):
            try:
                self.hash_table.search(pred)
                cessor_list.append(pred)
            except KeyError:
                pass
        return cessor_list

    # PAS METHODE PREDECESSEUR! CAR SI ENVOIE 010 et 101 ON EST FAITS!
    def _walk(self):
        copy = HashTabADN(self.hash_table)  # fait un mauvais resize mais cette table est temporaire anyway
        # TODO


class AlphabetError(Exception):
    pass


class ADNHashError(Exception):
    pass


class ADNCompressionError(Exception):
    pass


class Bucket:  # Creer qui si colision, sinon juste val direct (reduit besoin mem)
    col = 0

    @staticmethod
    def get_col():
        return Bucket.col

    @staticmethod
    def reset_col():
        Bucket.col = 0

    def __init__(self, val):
        self.list = []
        self.list.append(val)

    def add(self, val):
        self.list.append(val)
        Bucket.col += 1

    def __iter__(self):
        for ele in self.list:
            yield ele

    def __len__(self):
        return len(self.list)


class HashTabADN:

    def __init__(self, iterable=None, size=1000, word_length=21, factor=4):
        if iterable is None:
            print("Initial creation of table")

            self._size = size * factor  # length iterable * 2...
            self._used = 0
        else:
            print("Resized:")

            Bucket.reset_col()

            self._size = iterable.size() * factor
            self._used = len(iterable)

        print("\tsize: ", self._size)

        self._table = [None] * self._size

        self.load = self._used / self._size
        self._safe_bound = 4 ** word_length  # Nb total de mots dans notre langage

        self._prime = self._get_lowest_prime()
        self._scale = self._get_scale()

        if iterable is not None:
            for item in iterable:
                self.add(item, True)

    # Pour des questions d'optimisation du resize, on hash les nodes dans le BrujinGraph, a leur creation. De cette
    # facon, toutes les nodes sont deja hashees et on n'a qu'a les compresser pour les entrer dans la table, peu
    # importe si elles sont ajoutees pour la premiere fois a la table ou non
    def add(self, node, is_init=False):
        hashed = node.id

        com_hash = self._compress(hashed)
        old = self._table[com_hash]

        if isinstance(old, Node):
            bucket = Bucket(old)
            if hashed == old.id:
                return  # Ceci empeche le reisze si pas besoin ajouter...

            bucket.add(node)
        elif isinstance(old, Bucket):
            for item in old:  # On test si la node y est deja
                if item.hased == hashed:
                    return  # Ceci empeche le reisze si pas besoin ajouter...

            # Si non, on add
            old.add(node)
        else:
            self._table[com_hash] = node

        if not is_init:
            self._used += 1
            self._resize()

    def remove(self, name):
        node = self.search(name)

        node_to_refactor = None  # pour garder optimization pas bucket partout
        com_hash = self._hash_and_compress(name)
        temp_bucket = self._table[com_hash]
        if isinstance(temp_bucket, Bucket):
            temp_bucket.list.remove(node)

            if len(temp_bucket) == 1:
                node_to_refactor = temp_bucket.list[0]
                temp_bucket.list.remove(node_to_refactor)

        self._table[com_hash] = node_to_refactor  # None si pas besoin de refactor donc delete node, et sinon
        # remplace node :)

        self._used -= 1
        self._resize()

    def _resize(self):
        self.load = self._used / self._size
        if self.load >= 0.75:
            print("\tColisions: ", Bucket.get_col())
            resized = HashTabADN(self)
            self.__dict__.update(resized.__dict__)

    def __len__(self):
        return self._used

    def size(self):
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
        hashed = self.hash(name)  # prend hash non compresse
        item = self._table[self._compress(hashed)]  # on prend l'item correspondant dans la table

        if isinstance(item, Node) and item.id == hashed:  # Si est une Node, facile
            return item
        elif isinstance(item, Bucket):  # Si est Bucket,
            for ele in item:  # Cherche dans le bucket pour bonne clef
                if ele.id == hashed:
                    return ele

        raise KeyError

    """        for candidat in range(self._size, self._safe_bound):
            is_prime = True
            for facteur in range(2, candidat):
                if candidat % facteur == 0:
                    is_prime = False
                    break
            if is_prime:
                print("\tprime: ",  candidat)
                return candidat
        return max(self._size, self._safe_bound)"""

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
                print("\tprime: ", candidat)
                return candidat

    def _get_scale(self):
        low = self._prime // 2
        high = low + 1
        while high <= self._prime:
            if low % self._prime != 0:
                print("\tscale: ", low)
                return low
            if high % self._prime != 0:
                print("\tscale: ", high)
                return high

            low -= 1
            high += 1

        raise ADNCompressionError

    """
        key = ""
        for char in string:
            key += HashTabADN.conv_dict[char]

        try:
            key = int(key, 4)  # base 4! 4 letrres!
        except Exception:
            raise ADNHashError
    
        key = 0
        for index in range(len(string)):
            key += (_adn_to_int(string[-index]) * (10 ** index))
        try:
            key = int(str(key), 4)
        except Exception:
            raise ADNHashError
            
            
        try:
            key = int(str(self._r_builder(0, string)), 4)
        except Exception:
            raise ADNHashError       
            
    def _r_builder(self, key, string):
        if len(string) == 0:
            return key
        key = key*10 + self._adn_to_int(string[-1])
        return self._r_builder(key, string[:-1])
            """
    conv_dict = {'A': '0', 'C': '1', 'T': '2', 'G': '3'}

    def hash(self, string):
        key = ""
        for char in string:
            key += HashTabADN.conv_dict[char]

        try:
            key = int(key, 4)  # base 4! 4 letrres!
        except Exception:
            raise ADNHashError

        return key

    def _compress(self, key):
        if self._safe_bound > self._prime:
            key = key * self._scale  # M
            key = key % self._prime  # D
        return key % self._size  # Resize

    def _hash_and_compress(self, string):
        return self._compress(self.hash(string))

    def colision_test(self, str_tab):
        test_set = set()
        colisions = 0
        for string in str_tab:
            key = self._hash_and_compress(string)
            if key in test_set:
                colisions += 1
            else:
                test_set.add(key)

        print("Colisions:", colisions)
        return colisions
