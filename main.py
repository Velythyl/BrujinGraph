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