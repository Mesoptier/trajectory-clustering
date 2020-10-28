#!/usr/bin/python

import matplotlib.pyplot as plt

dataset = [ "Bladon Heath", "Church Hanborough", "Horspath", "Weston" ]

k = range(1,8)

l1 = [
(1.3697e+12, 1.3697e+12, 1.43467e+12),
(6.07749e+11, 6.07749e+11, 1.26408e+12),
(3.73933e+11, 3.73933e+11, 1.21088e+12),
(2.96759e+11, 2.96759e+11, 1.25193e+12),
(2.68214e+11, 2.67668e+11, 7.03523e+11),
(2.44655e+11, 2.38477e+11, 6.64298e+11),
(2.26044e+11, 2.13439e+11, 8.85437e+11)
]

l2 = [
(5.31548e+12, 5.31548e+12, 5.76173e+12),
(2.89599e+12, 2.89599e+12, 5.85325e+12),
(1.43642e+12, 1.43642e+12, 5.71193e+12),
(6.67904e+11, 6.67904e+11, 5.11298e+12),
(4.6676e+11, 4.6676e+11, 5.19374e+12),
(3.82112e+11, 3.82112e+11, 5.1595e+12),
(3.31524e+11, 3.31524e+11, 5.12202e+12)
]

l3 = [
(3.33577e+12, 3.33577e+12, 3.45371e+12),
(1.19382e+12, 1.19382e+12, 3.3683e+12),
(8.86647e+11, 8.86647e+11, 3.50996e+12),
(6.64148e+11, 7.36572e+11, 3.24332e+12),
(5.00386e+11, 5.31861e+11, 3.22719e+12),
(4.04839e+11, 4.04839e+11, 3.17699e+12),
(3.52352e+11, 3.52352e+11, 3.22048e+12)
]

l4 = [
(7.48491e+11, 7.48491e+11, 8.45139e+11),
(3.79356e+11, 3.79356e+11, 6.96003e+11),
(2.97515e+11, 2.97515e+11, 6.81738e+11),
(2.56122e+11, 2.71242e+11, 6.33506e+11),
(2.2729e+11, 2.2729e+11, 6.08575e+11),
(1.99449e+11, 2.12844e+11, 5.08472e+11),
(1.80787e+11, 1.86992e+11, 4.67532e+11)
]

data = [l1, l2, l3, l4]

labels = ["G.+PAM", "PAM", "G."]

def export_plot(X, Y_tuples, L, title, filename):
    font = {'family' : 'normal', 'size'   : 18}
    plt.rc('font', **font)
    #  rcParams.update({'figure.autolayout': True})

    fig,ax = plt.subplots()
    linestyles = ['-', '-.', '--']
    linewidths= [3,2,2]
    for i in range(3):
        Y = [ y[i] for y in Y_tuples ]
        objects = ax.plot(X, Y, linestyle=linestyles[i], alpha=1., label=L[i], linewidth=linewidths[i])

    plt.legend()
    plt.title(title)
    plt.xlabel("k")
    plt.ylabel("k-medians cost")
    plt.savefig(filename, bbox_inches='tight')
    plt.clf()

for i in range(4):
    export_plot(k, data[i], labels, dataset[i], dataset[i] + ".pdf")
