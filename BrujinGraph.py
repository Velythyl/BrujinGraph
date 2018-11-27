import gc
import math
import os
import sys
import time

import psutil


class Node:

    def __init__(self, string, id):
        self.kmer = string
        self.id = id


import multiprocessing as mp


class DeBrujinGraph:
    id_counter = 0

    @staticmethod
    def _build_id():
        temp = str(DeBrujinGraph.id_counter)
        while len(temp) < 10:
            temp = "" + temp

        DeBrujinGraph.id_counter += 1
        return temp

    def __init__(self, nodes, k=21, use_mp=False):
        self._alphabet = ['A', 'C', 'T', 'G']
        self._k = k
        self.hash_table = HashTabADN(size=int(len(nodes) * 8))
        # self.hash_table.colision_test(nodes)
        if use_mp is not False:
            self._mp_add_all(nodes)
        else:

            i = 0
            l = len(nodes)
            for name in nodes:
                for kmer in self._build_kmers(name):
                    self.add(kmer)
                print("\t", i, " of ", l, " added.")
                i += 1

    def _build_kmers(self, name):
        for i in range(len(name) - self._k + 1):
            yield name[i:i + self._k]

    def _split_in(self, list, n):
        split_list = []
        splitter = math.ceil(len(list) / n)
        for i in range(0, len(list), splitter):
            split_list.append(list[i:i + splitter])
        return split_list

    """     cpu_nb = mp.cpu_count()   # On ne veut pas les hyperthreads...
            nb_of_frags = cpu_nb
            split_list = self._split_in(nodes, nb_of_frags)
    
            pool = mp.Pool(cpu_nb)
            counter = 0
            for i in range(nb_of_frags):
                print("\tBuilding sub list", i + 1, "...")
                for node_list in pool.map(self._mp_hash_all_mapper, split_list[i]):
                    for node in node_list:
                        if node not in self.hash_table:
                            self.hash_table.add(node)
                    print("\t", counter, "of", len(split_list[i]), "added")
                    counter += 1
    
                pool.terminate()
                pool.join()
                gc.collect()
                pool = mp.Pool(cpu_nb)
    
                print("\t" + str(i + 1), "done of", nb_of_frags)
    
            pool.close()
            pool.join()
            
            
            cpu_nb = mp.cpu_count()   # On ne veut pas les hyperthreads...
        nb_of_frags = cpu_nb
        split_list = self._split_in(nodes, nb_of_frags)

        lock = mp.Lock()
        pool = mp.Pool(cpu_nb, initargs=(lock,))
        counter = 0
        for i in range(nb_of_frags):
            print("\tBuilding sub list", i + 1, "...")
            for node_list in pool.map(self._mp_hash_all_mapper, split_list[i]):
                for node in node_list:
                    if node not in self.hash_table:
                        self.hash_table.add(node)
                print("\t", counter, "of", len(split_list[i]), "added")
                counter += 1

            gc.collect()
            print("\t" + str(i + 1), "done of", nb_of_frags)

        pool.close()
        pool.join()
        
        cpu_nb = mp.cpu_count()   # On ne veut pas les hyperthreads...
        nb_of_frags = cpu_nb
        split_list = self._split_in(nodes, nb_of_frags)

        lock = mp.Lock()
        pool = mp.Pool(cpu_nb, initializer=self._mp_init, initargs=(lock,))
        memory_sentinels = [pool.apply(self._mp_get_sentinels) for i in range(cpu_nb)]
        this_pid = os.getpid()
        results = [pool.apply_async(self._mp_hash_all, args=(split_list[i],)) for i in range(nb_of_frags)]
        pool.close()
        counter = 0
        while True:
            total_mem_percent = psutil.Process(this_pid).memory_percent()
            for sentinel in memory_sentinels:
                total_mem_percent += psutil.Process(sentinel).memory_percent()

            if total_mem_percent > 20:
                lock.acquire()
                for res in results:
                    for node_list in res.get():
                        for node in node_list:
                            if node not in self.hash_table:
                                self.hash_table.add(node)
                        print("\t", counter, "of", len(nodes), "added")
                        counter += 1
                lock.release()

            if counter >= len(nodes):
                break

            time.sleep(1)

        pool.join()
        
        
        cpu_nb = mp.cpu_count()-1  # On ne veut pas les hyperthreads...

        lock = mp.Lock()
        pool = mp.Pool(cpu_nb, initializer=self._mp_init, initargs=(lock,))
        counter = 0

        # TODO LOCK ALL POOL WORKERS IF MEM >= X
        for node_list in pool.imap_unordered(self._mp_hash_all_mapper, nodes, chunksize=500):
            for node in node_list:
                if node not in self.hash_table:
                    self.hash_table.add(node)
            print("\t", counter, "of", len(nodes), "added")
            counter += 1

        pool.terminate()
        pool.join()
        gc.collect()
        pool = mp.Pool(cpu_nb)

        pool.close()
        pool.join()
            """

    import os
    import psutil
    def _mp_add_all(self, nodes):  # hash est op la plus couteuse, on la
        cpu_nb = mp.cpu_count() - 1  # On ne veut pas les hyperthreads...

        lock = mp.Lock()
        pool = mp.Pool(cpu_nb, initializer=self._mp_init, initargs=(lock,os.getpid(),))
        counter = 0

        iterator = pool.imap_unordered(self._mp_hash_all_mapper_test, nodes, chunksize=500)
        pool.close()
        for node_list in iterator:
            for node in node_list:
                if node not in self.hash_table:
                    self.hash_table.add(node)
            print("\t", counter, "of", len(nodes), "added")
            counter += 1

        pool.join()
        gc.collect()

    #https://stackoverflow.com/questions/28664720/how-to-create-global-lock-semaphore-with-multiprocessing-pool-in-python
    def _mp_init(self, lock_, main_pid_):
        global lock, main_pid
        lock = lock_
        main_pid = main_pid_

        """
        while True:
            
        for i in range(nb_of_frags):
            print("\tBuilding sub list", i + 1, "...")
            for node_list in pool.apply_async(self._mp_hash_all_mapper, split_list[i]):
                for node in node_list:
                    if node not in self.hash_table:
                        self.hash_table.add(node)
                print("\t", counter, "of", len(split_list[i]), "added")
                counter += 1

            gc.collect()
            print("\t" + str(i + 1), "done of", nb_of_frags)"""



    map_counter = 0
    def _mp_hash_all_mapper(self, node):
        kmer_list = []

        for kmer in self._build_kmers(node):
            kmer_list.append(Node(kmer, self._hash(kmer)))

        print("\tApprox.", DeBrujinGraph.map_counter*(mp.cpu_count()), "hashed")
        DeBrujinGraph.map_counter += 1

        return kmer_list

    def _mp_hash_all_mapper_test(self, node):
        kmer_list = []

        for kmer in self._build_kmers(node):
            kmer_list.append(Node(kmer, self._hash(kmer)))

        main_mem = psutil.Process(main_pid).memory_percent()

        print("\tApprox.", DeBrujinGraph.map_counter, "hashed, using", int(main_mem),"\b% memory in main process")
        DeBrujinGraph.map_counter += 1

        if main_mem >= 40:
            lock.acquire()
            while True:
                time.sleep(1)
                if psutil.Process(main_pid).memory_percent() <=25:  # TODO DEAD LOCK...
                    lock.release()
                    break



        return kmer_list

    def _mp_get_sentinels(self):
        time.sleep(1)   # S'assure d'avoir les bons pid
        return os.getpid()

    hash_counter = 0

    def _mp_hash_all(self, name_list):
        kmer_list = []

        for node in name_list:
            for kmer in self._build_kmers(node):
                kmer_list.append(Node(kmer, self._hash(kmer)))
            print("\tApprox.", DeBrujinGraph.hash_counter * (mp.cpu_count()), "hashed")
            DeBrujinGraph.hash_counter += 1
        return kmer_list

    def _hash(self, string):  # Sucre syntaxique
        return HashTabADN.hash(string)

    def __contains__(self, name: str) -> bool:  # détermine si le graphe de Brujin contient le noeud N
        try:
            self.hash_table.search_str(name)
            return True
        except KeyError:
            return False

    def __iter__(self):
        return self.nodes()  # retourne un itérable sur les noeuds du graphe TODO NOEUD OU STR???

    def load_factor(self) -> float:  # calcule le facteur de charge de la table de hachage sous-jacente
        return self.hash_table.load

    def add(self, kmer: str):  # ajoute le noeud N au graphe   TODO NOEUD OU STR
        node = Node(kmer, self._hash(kmer))
        if node not in self.hash_table:
            self.hash_table.add(node)

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
                self.hash_table.search_str(pred)
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

    def __init__(self, size=400000, word_length=21):
        print("\tInitial creation of table:")

        self._size = size  # length iterable * 2...
        self._used = 0

        print("\t\t\tsize: ", self._size)

        self._table = [None] * self._size

        self.load = self._used / self._size
        self._safe_bound = 4 ** word_length  # Nb total de mots dans notre langage

        self._prime = 92821  # self._get_lowest_prime()
        self._scale = 46410 # self._get_scale()

    def _rebuild(self, old, old_used, factor=4):
        print("\t\tResized: rebuilding...")

        Bucket.reset_col()

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
        hashed = node.id

        com_hash = self._compress(hashed)
        old = self._table[com_hash]

        if isinstance(old, Node):
            bucket = Bucket(old)
            bucket.add(node)
            self._table[com_hash] = bucket
        elif isinstance(old, Bucket):
            old.add(node)
        else:
            self._table[com_hash] = node

        if not is_init:
            self._used += 1
            self._resize()

    def remove(self, name):
        node = self.search_str(name)

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
            print("\t\t\tColisions: ", Bucket.get_col())
            self._rebuild(self, self._used)
            """
            resized = HashTabADN(self)
            self.__dict__.update(resized.__dict__)"""

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
            if ele is None:
                continue
            elif isinstance(ele, Node):
                yield ele
            elif isinstance(ele, Bucket):
                for sub_ele in ele:
                    yield sub_ele
            else:
                raise Exception

    def search_node(self, node):
        item = self._table[self._compress(node.id)]  # on prend l'item correspondant dans la table

        if isinstance(item, Node) and item.id == node.id:  # Si est une Node, facile
            return item
        elif isinstance(item, Bucket):  # Si est Bucket,
            for ele in item:  # Cherche dans le bucket pour bonne clef
                if ele.id == node.id:
                    return ele

        raise KeyError

    def search_str(self, name) -> Node:
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

    @staticmethod
    def hash(string):
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

        print("\tColisions:", colisions)
        return colisions
