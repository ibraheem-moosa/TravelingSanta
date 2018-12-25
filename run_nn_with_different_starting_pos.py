import random
import os
N = 197770
for i in range(10):
    sp = random.randint(0, N)
    os.system("./a.out cities.csv result-{}.csv {}".format(sp, sp))
