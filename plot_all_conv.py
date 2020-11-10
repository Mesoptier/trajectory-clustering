from plot_clustering import *

for i in range(10):
    plot_clustering
    filename = "convergence/Bladon & Church route recapping/bladon heath/brc_" + str(2 * i) + "_5"
    figure_name = "brc_" + str(2 * i) + "_5"
    radius = .02
    file = open(filename, "r")
    plot_clustering(file, figure_name, radius)
