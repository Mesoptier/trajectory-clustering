import matplotlib.pyplot as plt
import pylab as pl
import numpy as np
import os
import pandas
import sys
from matplotlib import collections  as mc


path_file = open('warping_path.txt').read().rstrip().split("\n")
c1 = open('c1.txt').read().rstrip().split("\n")
c2 = open('c2.txt').read().rstrip().split("\n")

c1_dist = []
c2_dist = []

c1_total = 0
c2_total = 0

y_min = 1000000
y_max = -1000000
x_min = 1000000
x_max = -1000000

for i in range(1, len(c1)):
    p1 = c1[i].split(",")
    p2 = c1[i-1].split(",")
    dist = pow(
        pow(float(p1[0])-float(p2[0]), 2) + pow(float(p1[1])-float(p2[1]), 2)
    , 0.5)
    c1_dist.append(c1_total + dist)
    print(c1_dist[len(c1_dist) - 1])
    c1_total += dist

    if float(p1[0]) > x_max:
        x_max = float(p1[0])
    if float(p1[0]) < x_min:
        x_min = float(p1[0])

    if float(p1[1]) > y_max:
        y_max = float(p1[1])
    if float(p1[1]) < y_min:
        y_min = float(p1[1])


    if float(p2[0]) > x_max:
        x_max = float(p2[0])
    if float(p2[0]) < x_min:
        x_min = float(p2[0])

    if float(p2[1]) > y_max:
        y_max = float(p2[1])
    if float(p2[1]) < y_min:
        y_min = float(p2[1])


for i in range(1, len(c2)):
    p1 = c2[i].split(",")
    p2 = c2[i-1].split(",")
    dist = pow(
        pow(float(p1[0])-float(p2[0]), 2) + pow(float(p1[1])-float(p2[1]), 2)
    , 0.5)
    c2_dist.append(c2_total + dist)
    print(c2_dist[len(c2_dist) - 1])
    c2_total += dist

    if float(p1[0]) > x_max:
        x_max = float(p1[0])
    if float(p1[0]) < x_min:
        x_min = float(p1[0])

    if float(p1[1]) > y_max:
        y_max = float(p1[1])
    if float(p1[1]) < y_min:
        y_min = float(p1[1])


    if float(p2[0]) > x_max:
        x_max = float(p2[0])
    if float(p2[0]) < x_min:
        x_min = float(p2[0])

    if float(p2[1]) > y_max:
        y_max = float(p2[1])
    if float(p2[1]) < y_min:
        y_min = float(p2[1])


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

seg_c = mc.LineCollection(path_segs, linewidths=2, color='w')

lc = mc.LineCollection(data, linewidths=2, linestyle='dashed', color = "g")
fig, (ax1, ax2) = plt.subplots(1,2, figsize=(20,8))
ax1.add_collection(lc)
ax1.add_collection(seg_c)
ax1.autoscale(False)
ax1.set_xlim(0, c2_total)
ax1.set_ylim(0, c1_total)
ax1.margins(0.1)

heatmap_data = np.loadtxt("heat_map.dat")
# ax1.imshow(heat_map, cmap='hot', interpolation='nearest')

min_value = heatmap_data.min()
max_value = heatmap_data.max()

ex = [0 , c2_total, 0 , c1_total]

heatmap = ax1.imshow(heatmap_data, cmap='hot', vmin=min_value, vmax=max_value, origin='lower', extent=ex)
ax1.set_aspect(c1_total / c2_total)
# cbar = plt.colorbar(heatmap, shrink=1., pad=0.)
# cbar.set_label(r'Fr\'echet Distance' , fontsize=40)
# plt.savefig(output_filename, bbox_inches='tight')

# color axes

y_ax = [[
    (0, 0), (0, c1_total)
]]

x_ax = [[
    (0, 0), (c2_total, 0)
]]

y_c = mc.LineCollection(y_ax, linewidths=8, color='y')
x_c = mc.LineCollection(x_ax, linewidths=8, color='m')

ax1.add_collection(y_c)
ax1.add_collection(x_c)


# plotting curves in image space

edges1 = []
edges2 = []

for i in range(1, len(c1)):
    p = c1[i-1].split(",")
    q = c1[i].split(",")
    edges1.append(
        [(float(p[0]), float(p[1])), (float(q[0]), float(q[1]))]
    )


for i in range(1, len(c2)):
    p = c2[i-1].split(",")
    q = c2[i].split(",")
    edges2.append(
        [(float(p[0]), float(p[1])), (float(q[0]), float(q[1]))]
    )

im_c1 = mc.LineCollection(edges1, linewidths=2, color="y")
im_c2 = mc.LineCollection(edges2, linewidths=2, color="m")
ax2.add_collection(im_c1)
ax2.add_collection(im_c2)
ax2.autoscale(False)
ax2.set_xlim(x_min, x_max)
ax2.set_ylim(y_min, y_max)
ax2.margins(0.1)

ax2.set_aspect(1)

c1_vert = {
    "x": [],
    "y": []
}

for e in c1:
    e = e.split(",")
    c1_vert["x"].append(float(e[0]))
    c1_vert["y"].append(float(e[1]))

c1_vert = pandas.DataFrame(c1_vert)

c1_vert.plot(kind="scatter", x="x", y="y", color="black", ax = ax2)

c2_vert = {
    "x": [],
    "y": []
}

for e in c2:
    e = e.split(",")
    c2_vert["x"].append(float(e[0]))
    c2_vert["y"].append(float(e[1]))

c2_vert = pandas.DataFrame(c2_vert)

c2_vert.plot(kind="scatter", x="x", y="y", color="black", ax = ax2)


print([x_min, x_max, y_min, y_max])

fig.savefig('warping_path.png')