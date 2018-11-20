# test de colisions avec notre fonction de hashage
import itertools

alphabet = ['A', 'T', 'C', 'G']
str_tab = []
breaker = 0
for i in itertools.product(alphabet, repeat = 6):


    str_tab.append(''.join(i))

from HashTab import HashTabADNSeed2
temp = HashTabADNSeed2()
temp.colision_test(str_tab)
