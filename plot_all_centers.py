#!/usr/bin/env python3

import pathlib
import matplotlib.pyplot as plt

sites = [['a55', 'brc', 'c17', 'c35', 'p29', 'p39', 'p94'],
    ['a94', 'c22', 'c70', 'k77', 'l29', 'liv', 'r47', 's93'],
    ['H22', 'H27', 'H30', 'H35', 'H38', 'H41', 'H42', 'H71'],
    ['H23', 'H31', 'H32', 'H34', 'H36', 'H50', 'H58', 'H62']]
paths = ['Bladon & Church route recapping/bladon heath',
    'Bladon & Church route recapping/church hanborough',
    'Horspath', 'Weston']
names = ['bladon', 'church', 'horspath', 'weston']
centers = []
base = pathlib.Path('convergence')

for (i, p) in enumerate(paths):
    centers.append([])
    for s in sites[i]:
        all_pigeons_d = sorted((base / p).rglob(s + '_?_5'))
        if not all_pigeons_d:
            continue
        with open(all_pigeons_d[-2], 'r') as file:
            file.readline()
            number_of_curves = int(file.readline())
            for j in range(number_of_curves):
                file.readline()
            centers[i].append(list(map(float,
                file.readline().rstrip().split(' '))))

for (i, cs) in enumerate(centers):
    fig = plt.figure(dpi=600, frameon=False)
    ax = fig.add_subplot(111)
    for center in cs:
        for j in range(0, len(center) - 2, 2):
            plt.plot((center[j], center[j + 2]), (center[j + 1], center[j + 3]),
                'm', linewidth=.7)
            ax.add_artist(plt.Circle((center[j], center[j + 1]), .02,
                color='m'))
        ax.add_artist(plt.Circle((center[-2], center[-1]), .02, color='m'))
    ax.set_aspect('equal', adjustable='box')
    plt.axis('off')
    plt.savefig('figures/' + names[i] + '_centers.png', dpi=600,
        bbox_inches='tight')
