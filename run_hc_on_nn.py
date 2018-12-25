import os
import sys

tours = os.listdir(sys.argv[1])
for t in tours:
    print(t)
    os.system("./a.out cities.csv {} {}".format(os.path.join(sys.argv[2], t), os.path.join(sys.argv[1], t)))
