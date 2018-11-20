import random

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

class HashTabADN:

    def __init__(self):
        table = []
        self.size = 5120
        load = 0.0

    def _hash_it(self, string):
        key = 0
        for index in range(len(string)):
            key += (_adn_to_int(string[index]) * (10 ** index))
        return key

    def _compress(self, key):
        return key % self.size

    def hash(self, string):
        return self._compress(self._hash_it(string))

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
