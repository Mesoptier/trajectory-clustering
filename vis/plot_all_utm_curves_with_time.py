#!/usr/bin/python3
# Gets a directory as argument and plots all the curve files that are in there.

import matplotlib
import matplotlib.pyplot as plt
import sys
import glob

prefix = "../data/with_utm/Data_for_Mann_et_al_RSBL/"
utm_dirs = [ "Bladon & Church route recapping/church hanborough/c70/utm","Bladon & Church route recapping/church hanborough/r47/utm","Bladon & Church route recapping/church hanborough/s93/utm","Bladon & Church route recapping/church hanborough/liv/utm","Bladon & Church route recapping/church hanborough/a94/utm","Bladon & Church route recapping/church hanborough/l29/utm","Bladon & Church route recapping/church hanborough/c22/utm","Bladon & Church route recapping/church hanborough/k77/utm","Bladon & Church route recapping/bladon heath/c35/utm","Bladon & Church route recapping/bladon heath/p29/utm","Bladon & Church route recapping/bladon heath/c17/utm","Bladon & Church route recapping/bladon heath/brc/utm","Bladon & Church route recapping/bladon heath/a55/utm","Bladon & Church route recapping/bladon heath/p94/utm","Bladon & Church route recapping/bladon heath/p39/utm","Weston/H58/utm","Weston/H23/utm","Weston/H36/utm","Weston/H32/utm","Weston/H34/utm","Weston/H31/utm","Weston/H50/utm","Weston/H62/utm","Horspath/H42/utm","Horspath/H22/utm","Horspath/H30/utm","Horspath/H41/utm","Horspath/H71/utm","Horspath/H38/utm","Horspath/H27/utm","Horspath/H35/utm" ]

def readCurves(curve_directory):
    curves = []

    for filename in sorted(glob.glob(curve_directory + "/*")):
        curve = [[],[]] # a curve is an x and a y list
        if filename == curve_directory + "/dataset.txt":
            continue
        with open(filename, 'r') as f:
            f.readline() # read header
            for line in f:
                A = line.split()
                if A:
                    curve[0].append(float(A[0]))
                    curve[1].append(float(A[1]))
        curves.append(curve)

    return curves

def saveCurves(curves, filename):
    values = [ (1./len(curves))*i for i in range(len(curves)) ]
    for curve,value in zip(curves,values):
        plt.plot(curve[0], curve[1], c=(value,0.,0.))
    plt.axis('off')
    plt.savefig(filename, bbox_inches='tight')
    plt.clf()

if len(sys.argv) != 1:
    print("USAGE: ./plot_curves.py")
    quit()

for i,utm_dir in enumerate(utm_dirs):
    curve_directory = prefix + utm_dir
    filename = "out/" + str(i) + ".pdf"
    curves = readCurves(curve_directory)
    saveCurves(curves, filename)
