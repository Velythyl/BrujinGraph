from HashTab import HashTabADN


class Node:
    def __init__(self, string, hash):
        self.name = string
        self.hash = hash
        self.predecessors = set()
        self.successors = set()


class DeBrujinGraph:

    def __init__(self, names: Iterable[str], k=21):
        self.hash_table = HashTabADN()
        self._alphabet = ['A', 'C', 'T', 'G']
        for name in names:
            self.add(name)

    def _hash(self, string):  # Sucre syntaxique
        return self.hash_table.hash(string)

    def __contains__(self, name: str) -> bool:  # détermine si le graphe de Brujin contient le noeud N
        try:
            self.hash_table.search(name)
            return True
        except KeyError:
            return False

    def __iter__(self) -> Iterable[str]:
        return self.nodes()  # retourne un itérable sur les noeuds du graphe TODO NOEUD OU STR???

    def load_factor(self) -> float:  # calcule le facteur de charge de la table de hachage sous-jacente
        return self.hash_table.load

    def add(self, name: str):  # ajoute le noeud N au graphe   TODO NOEUD OU STR

        node = Node(name, self.hash_table.hash(name))
        self.hash_table.add(node)

        for candidate in self._get_pred(name):
            try:
                pred = self.hash_table.search(candidate)
                self._link(pred, node)
            except Exception:
                pass

        for candidate in self._get_succ(name):
            try:
                succ = self.hash_table.search(candidate)
                self._link(node, succ)
            except Exception:
                pass

    def _link(self, pred, succ):
        pred.successors.add(succ)
        succ.predecessors.add(pred)

    def _get_succ(self, name: str):
        name = name[1:]
        for char in self._alphabet:
            yield name + char

    def _get_pred(self, name: str):
        name = name[:-1]
        for char in self._alphabet:
            yield char + name

    def remove(self, name: str):  # enlève le noeud N du graphe
        self.hash_table.remove(name)

    def nodes(self) -> Iterable[str]:  # retourne un itérable sur les noeuds du graphe TODO NOEUD OU STR???
        for node in self.hash_table:
            yield node.name

    def predecessors(self, name: str) -> Iterable[str]:  # retourne tous les prédécesseur du noeud N
        return self.hash_table.search(name).predecessors

    def successors(self, name: str) -> Iterable[str]:  # retourne tous les successeurs du noeud N
        return self.hash_table.search(name).successors

    # PAS METHODE PREDECESSEUR! CAR SI ENVOIE 010 et 101 ON EST FAITS!
    def _walk(self):
        copy = HashTabADN(self.hash_table) # fait un mauvais resize mais cette table est temporaire anyway
            #TODO


