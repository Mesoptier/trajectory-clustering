#!/usr/bin/python3
# Gets a directory as argument and plots all the curve files that are in there.

import matplotlib
import matplotlib.pyplot as plt
import sys
import glob

def readCurves(curve_directory):
    curves = []

    for filename in glob.glob(curve_directory + "/*"):
        curve = [[],[]] # a curve is an x and a y list
        if filename == curve_directory + "/dataset.txt":
            continue
        with open(filename, 'r') as f:
            f.readline() # read header
            for line in f:
                A = line.split()
                if A:
                    curve[0].append(A[0])
                    curve[1].append(A[1])
        curves.append(curve)

    return curves

def showCurves(curves):
    for curve in curves:
        plt.plot(curve[0], curve[1])
    plt.show()

if len(sys.argv) != 2:
    print("USAGE: ./plot_curves.py <curve_directory>")
    quit()

curve_directory = sys.argv[1]
curves = readCurves(curve_directory)
showCurves(curves)
