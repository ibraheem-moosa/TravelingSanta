import sys
from matplotlib import collections as mc
import pylab as pl

with open(sys.argv[1]) as tour_file, open(sys.argv[2]) as cities_file:
    tour = [l for l in tour_file]
    # ignore first
    tour = tour[1:]
    tour = list(map(int, tour))

    cities = [l for l in cities_file]
    # ignore first line
    cities = cities[1:]
    cities = [l.split(',') for l in cities]
    cities = dict([(int(l[0]), (float(l[1]), float(l[2]))) for l in cities])

    lines = [[cities[tour[i]], cities[tour[i+1]]] for i in range(len(tour) - 1)]
    lc = mc.LineCollection(lines, linewidth=2)
    fig, ax = pl.subplots(figsize=(20,20))
    ax.set_aspect('equal')
    ax.add_collection(lc)
    ax.autoscale()
    pl.show()

