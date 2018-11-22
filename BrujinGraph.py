from HashTab import HashTabADN

class Node:
    def __init__(self, string):
        self.name = string
        self.hash = HashTabADN.hash(string)
        self.predecessors = []
        self.successors = []

class DeBrujinGraph:




    def __init__(self, names: Iterable[str], k=21):
        self.hash_table = HashTabADN()
        for name in names:
            self.add(name)

    def _hash(self, string):    # Sucre syntaxique
        return self.hash_table.hash(string)

    def __contains__(self, name: str) -> bool: # détermine si le graphe de Brujin contient le noeud N
        try:
            self.hash_table.search(name)
            return True
        except KeyError:
            return False

    def __iter__(self) -> Iterable[str]:
        return self.nodes() # retourne un itérable sur les noeuds du graphe TODO NOEUD OU STR???

    def load_factor(self) -> float: # calcule le facteur de charge de la table de hachage sous-jacente
        return self.hash_table.load

    def add(self, name: str):   # ajoute le noeud N au graphe   TODO NOEUD OU STR
        # TODO creer node selon str...
        self.hash_table.add(node)

    def remove(self, name: str):    # enlève le noeud N du graphe
        self.hash_table.remove(name)

    def nodes(self) -> Iterable[str]:   # retourne un itérable sur les noeuds du graphe TODO NOEUD OU STR???
        for node in self.hash_table:
            yield node.name

    def predecessors(self, name: str) -> Iterable[str]:     # retourne tous les prédécesseur du noeud N
        return self.hash_table.search(name).predecessors

    def successors(self, name: str) -> Iterable[str]:       # retourne tous les successeurs du noeud N
        return self.hash_table.search(name).successors