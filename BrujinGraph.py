import pickle
import sys

sys.setrecursionlimit(10000)


# Opere sur des noms/kmer/str

class DeBrujinGraph:
    def __init__(self, nodes=None, k=21):
        self._k = k
        self.hash_table = HashTabADN()
        if nodes is not None:  # Si on cree un graph avec un iterateur
            for kmer in nodes:
                self.add(kmer)

    # détermine si le graphe de Brujin contient le noeud N
    def __contains__(self, N: str) -> bool:
        return N in self.hash_table  # On test simplement si la table contient N

    # retourne un itérable sur les noeuds du graphe
    def __iter__(self):
        return self.nodes()

    # On utilise pas cette fonction puisqu'on a code une table de hashage sous-jacente, mais elle etait demandee...
    def load_factor(self) -> float:  # calcule le facteur de charge de la table de hachage sous-jacente
        return self.hash_table.load

    # ajoute le noeud N au graphe
    def add(self, N: str):
        self.hash_table.add(hash(N))  # On ajoute simplement le hash du noeud a la table

    # ajoute les noeud pre-hashes de node_list au graphe
    def add_hashed_node_list(self, node_list):
        for node in node_list:
            self.hash_table.add(node)

    # Retourne les successeurs d'un kmer (ceci ne teste pas si ils appartiennent au graph)
    def _get_succ(self, kmer: str):
        name = kmer[1:]
        for char in ['A', 'T', 'C', 'G']:
            yield name + char

    # Retourne les predecesseurs d'un kmer (ceci ne teste pas si ils appartiennent au graph)
    def _get_pred(self, kmer: str):
        name = kmer[:-1]
        for char in ['A', 'T', 'C', 'G']:
            yield char + name

    # enlève le noeud N du graphe
    def remove(self, N: str):
        self.hash_table.remove(hash(N))  # On enleve le hash du noeud de la table

    # retourne un itérable sur les noeuds du graphe
    def nodes(self):
        for node in self.hash_table:  # Pour chaque node dans la table
            yield unhash(node)  # On la de-hash et on la yield

    # retourne tous les prédécesseur connus du noeud N
    def predecessors(self, N: str):
        return self.cessors(N, True)

    # retourne tous les successeurs connus du noeud N
    def successors(self, N: str):
        return self.cessors(N, False)

    # Retourne les succ ou prede -cesseurs connus, selon le bool is_pred
    def cessors(self, kmer: str, is_pred: bool):
        cessor_list = []
        for item in self._get_pred(kmer) if is_pred else self._get_succ(kmer):  # Pour chaque cesseur possible
            if item in self.hash_table:  # On regarde s'il est dans la table
                cessor_list.append(item)  # Si oui, on l'ajoute a la liste de cesseurs
        return cessor_list  # Puis, on retourne la liste

    # Iterateur qui retourne tous les departs pour un parcours de graph
    def get_all_starts(self):
        for ele in self.hash_table:  # Pour chaque noeud
            temp = unhash(ele)
            if len(self.predecessors(temp)) == 0:  # On teste s'il a des predecesseurs
                yield temp  # Si non, on le yield

    # Parcours le graph en entier d'apres ses starts, les noeuds sans parents
    def walk(self):
        for start in self.get_all_starts():  # Pour chaque start
            temp_set = set(start)  # On initialise un set
            for contig in self._walk(start, temp_set, start):  # On marche depuis ce start
                yield contig  # On yield les contigs resultant de cette marche

    # Parcours le graph a partir d'un noeud at, "vers le bas" (seulement vers les successeurs)
    def _walk(self, at, closed, contig):
        succ_list = self.successors(at)  # On prend les successeurs connus de at
        if len(succ_list) == 0:  # S'il n'y en a pas,
            yield contig  # Parcours finit
        else:  # Sinon
            for succ in self.successors(at):  # Pour chaque successeur
                if succ in closed:
                    yield contig + succ[-1]  # S'il a deja ete visite (loop), on yield le contig approprie
                else:
                    closed.add(succ)  # Sinon, on le marque comme visite
                    yield from self._walk(succ, closed, contig + succ[-1])  # Et on yield du prochain appel recursif

    def save(self, f="DBG.gra"):
        with open(f, "wb") as file:
            pickle.dump(self.__dict__, file)

    @staticmethod
    def loader(f="DBG.gra"):
        with open(f, "rb") as file:
            return pickle.load(file)


# Iterateur qui retourne tous les kmers pouvant etre construits d'apres un string "name"
def build_kmers(name, k=21):
    for i in range(len(name) - k + 1):  # Pour chaque index du string et tant qu'on a k-1 char a droite
        yield name[i:i + k]  # On retourne l'index + les k-1 prochains char


# Opere sur des int/hash/compressions

class HashTabADN:

    def __init__(self, size=400000):
        print("\tInitial creation of table:")

        self._size = size  # Au debut, on a une grosseur arbitraire
        self._used = 0  # Et 0 places dans la table d'utilisees

        print("\t\t\tsize: ", self._size)

        self._table = [None] * self._size  # On cree une table vide de la grosseur demandee

        self.load = 0.0  # On utilise rien au debut

        self._prime = 92821  # Prime pour MAD
        self._scale = 46410  # scale pour MAD

    def _rebuild(self, old, old_used, factor=4):
        print("\t\tResized: rebuilding...")

        self._size = old_used * factor  # On augmente la table par un facteur donne
        print("\t\t\tNew size:", self._size)
        self._used = old_used  # On utilise le meme nombre d'emplacements qu'avant
        self._table = [None] * self._size  # On reconstruit la table a la bonne grandeur
        self.load = self._used / self._size  # On a un load connu

        for item in old:  # Pour chaque item de la table d'avant
            self.add(item, True)  # On l'ajoute a celle-ci

    # Ajoute une node deja hashee a la table
    #
    # Si rebuilding est False, on est dans un cas normal et on gere le used et le load
    # Sinon, on est en train de rebuild et donc on a pas besoin de le faire
    #
    # Pour des questions d'optimisation du resize, on hash les nodes dans le BrujinGraph, a leur creation. De cette
    # facon, toutes les nodes sont deja hashees et on n'a qu'a les compresser pour les entrer dans la table, peu
    # importe si elles sont ajoutees pour la premiere fois a la table ou non
    def add(self, node, rebuilding=False):
        if node in self:  # Si la node est deja dans la table
            return  # On s'arrete

        com_hash = self._compress(node)  # Sinon, on prend la compression de la node
        old = self._table[com_hash]  # Et on va voir ce qu'il y a a cet emplacement dans la table

        if old is None:  # S'il n'y a rien, on peut simplement y assigner la node
            self._table[com_hash] = node
        elif isinstance(old, list):  # S'il y a une liste (un bucket), on y rajoute la node
            old.append(node)
        else:  # Sinon, on y a deja une node. On la remplace par un bucket de cette
            self._table[com_hash] = [old, node]  # node et de celle qu'on est en train d'ajouter

        if not rebuilding:  # Si on est pas en train de rebuild
            self._used += 1  # On a un used de plus
            self.load = self._used / self._size  # Et le load change
            if self.load >= 0.75:
                self._rebuild(self, self._used)  # Si le load grimpe >= a 0.75, on rebuild la table

    # Enleve le hash d'une node de la table.
    def remove(self, node):
        if node not in self:  # Si la node n'est pas dans la table, on s'arrete
            return

        self._used -= 1  # Sinon, on a un used de moins
        self.load = self._used / self._size  # Et le load change sans chance de resize

        node_to_refactor = None  # node_to_refactor est assume None
        com_hash = self._compress(node)  # On va chercher la compression de la node

        old = self._table[com_hash]  # Et ce qu'il y avait a cette compression
        if isinstance(old, list):  # Si c'etait une liste
            old.remove(node)  # On y enleve la node

            if len(old) == 1:  # Pour garder l'optimisation d'avoir des listes seulement si necessaire,
                node_to_refactor = old[0]  # si la liste n'a qu'un item, on met cet item dans node_to_refactor
            else:
                return  # Sinon, finit!

        self._table[com_hash] = node_to_refactor  # Si on avait pas de liste, node_to_refactor est None et on enleve
        # simplement la node de la table. Mais sinon, on la remplace en fait par l'item qui se retrouvait seul dans la
        # liste!

    # Retourne la grandeur reelle de la table, son size
    def size(self):
        return self._size

    # Retourne le nombre d'elements dans la table, son used
    def __len__(self):
        return self._used

    # Indique si node se trouve dans la table
    def __contains__(self, node):
        try:
            if isinstance(node, str):  # Parfois il est utilie de pouvoir demander si un string se trouve dans
                self.search_node(hash(node))  # la table depuis le graph. Donc, on hash simplement le string puis on
            else:  # cherche ce hash.
                self.search_node(node)  # Sinon, on cherche simplement le hash
            return True  # Si la recherche reussit, on se rend ici et on retourne True
        except KeyError:
            return False  # Sinon, le try-except nous amene ici et on retourne False

    # Iterateur qui yield tous les elements de la table
    def __iter__(self):
        for ele in self._table:  # Pour chaque ele directement dans la table
            try:  # On essaie de retourne tous les elements de ele comme si ele etait une liste
                for ele2 in ele:
                    yield ele2
            except Exception:  # Si ele n'etait pas une liste,
                if ele is None:
                    continue
                yield ele  # On retourne ele lui-meme s'il n'est pas None

    # Cherche la table pour une node, et retourne la node si on la trouve
    def search_node(self, node) -> int:
        com_hash = self._compress(node)  # On compresse pour trouver l'emplacement
        item = self._table[com_hash]  # On prend l'item correspondant dans la table

        if item == node:  # Si est une Node et que c'est celle qu'on cherche,
            return com_hash  # on la retourne

        try:  # Sinon, on suppose que item est une liste
            if node in item:
                return node  # On retourne node si node est dans cette liste
        except Exception:
            pass

        raise KeyError  # Si tout echoue, on lance un KeyError

    # Compresse une key pour qu'elle "rentre" dans la table
    # On ne fait pas MAD complet: en effet, les mots possibles sont tellement nombreux qu'on s'attend a avoir une tres
    # bonne distribution rien qu'avec D. De plus, prendre des primes appropries pour la grosseur de la table serait un
    # probeleme assez couteux...
    def _compress(self, key):
        return key % self._size  # Resize


conv_dict = {'A': '0', 'C': '1', 'T': '2', 'G': '3', '0': 'A', '1': 'C', '2': 'T', '3': 'G'}


# Retourne le string duquel key provient (defait un hashage)
def unhash(key):
    # Passe de la base 10 a la base 4
    # loosely inspired de https://www.codevscolor.com/python-convert-decimal-ternarybase-3/
    def _to_quaternary(num):
        q, r = divmod(num, 4)  # On divise le nombre par 4, et on garde le reste aussi.
        if q == 0:  # 4         # S'il n'y a pas de q (donc r rentre dans 4)
            return str(r)  # on retourne r
        else:  # Sinon,
            return _to_quaternary(q) + str(r)  # on retourne la base 4 de q, + r

    key = _to_quaternary(key)  # On prend la base 4 de la clef
    string = ""
    for char in key:  # Pour chaque char numerique (0,1,2,3) dans la clef en base 4,
        string += conv_dict[char]  # On ajoute la lettre correspondante (voir conv_dict) a la fin d'un string

    while len(string) < 21:  # Si la clef commencait par des A, le debut de la clef manque des 0! Donc, on ajoute
        string = 'A' + string  # A, la lettre correspondant a 0, jusqu'a ce que la clef mesure 21 de long.

    # On ne teste pas pour voir si _to_quaternary plante ici: en effet, comme toutes les nodes devraient etre passees
    # par hash AVANT de se faire unhash, on est certains que ca va marcher

    return string


class ADNHashError(Exception):  # Exception lancee lorsque hash ne fonctionne pas (probablement parce qu'on lui
    pass  # demande de hasher quelque chose n'etant pas compose que par A, T, C, G


# Hashe un string compose des lettres A, T, C, G
def hash(string):
    try:
        key = ""
        for char in string:  # Pour chaque char dans le string,
            key += conv_dict[char]  # on ajoute sa correspondance en chiffre (voir conv_dict) a un string "key"

        # Comme notre alphabet de comprends que 4 lettres, le string converti en la correspondance numerique de ces
        # lettres est en fait un string de numbre en base 4!
        key = int(key, 4)  # Donc, on transforme la clef en un int en base 4 (ceci sauve de l'espace en memoire)

        return key
    except Exception:
        raise ADNHashError  # Si quelque chose plante, on retourne une exception personnalisee pour l'indiquer.
