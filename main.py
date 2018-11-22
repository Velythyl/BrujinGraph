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
import pickle
upper_bound = 4**21
prime_list = [2, 3]
for candidat in range(4, upper_bound):
    is_prime = True
    for facteur in range(2, candidat):
        if candidat % facteur == 0:
            is_prime = False
            break
    if is_prime:
        print(candidat)
        prime_list.append(candidat)

with open('prime_list.pkl', 'wb') as fp:
    pickle.dump(prime_list, fp)