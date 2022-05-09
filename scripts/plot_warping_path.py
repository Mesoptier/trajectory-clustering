import matplotlib as plt
import pylab as pl
import os
import sys
from matplotlib import collections  as mc

path_file = open('warping_path.txt').read().rstrip().split("\n")
c1 = open('c1.txt').read().rstrip().split("\n")
c2 = open('c2.txt').read().rstrip().split("\n")

c1_dist = []
c2_dist = []

c1_total = 0
c2_total = 0

for i in range(1, len(c1)):
    p1 = c1[i].split(",")
    p2 = c1[i-1].split(",")
    dist = pow(
        pow(float(p1[0])-float(p2[0]), 2) + pow(float(p1[1])-float(p2[1]), 2)
    , 0.5)
    c1_dist.append(c1_total + dist)
    print(c1_dist[len(c1_dist) - 1])
    c1_total += dist

for i in range(1, len(c2)):
    p1 = c2[i].split(",")
    p2 = c2[i-1].split(",")
    dist = pow(
        pow(float(p1[0])-float(p2[0]), 2) + pow(float(p1[1])-float(p2[1]), 2)
    , 0.5)
    c2_dist.append(c2_total + dist)
    print(c2_dist[len(c2_dist) - 1])
    c2_total += dist


data = []

for i in range(len(c1_dist)):
    data.append(
        [(0, c1_dist[i]), (c2_total, c1_dist[i])]
    )

for i in range(len(c2_dist)):
    data.append(
        [(c2_dist[i], 0), (c2_dist[i], c1_total)]
    )

print(data)

data.append([(0, 0), (0, c1_total)])
data.append([(0, c1_total), (c2_total, c1_total)])
data.append([(0, 0), (c2_total, 0)])
data.append([(c2_total, 0), (c2_total, c1_total)])

path_segs = []

for line in path_file:
    coords = line.split(" ")
    path_segs.append([
        (float(coords[0]), float(coords[1])), (float(coords[2]), float(coords[3]))
    ])

seg_c = mc.LineCollection(path_segs, linewidths=2, color='r')

lc = mc.LineCollection(data, linewidths=1, linestyle='dashed')
fig, ax = pl.subplots()
ax.add_collection(lc)
ax.add_collection(seg_c)
ax.autoscale()
ax.margins(0.1)
fig.savefig('warping_path.png')