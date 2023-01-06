import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def plot_results(dataset):
    results = pd.read_csv('2d_paper_results/running_time_results_' + dataset + '.csv')
    pieces = results["number_of_pieces"].to_list()
    time = results["time"].to_list()

    x = results["i"].to_list()

    plt.plot(x, pieces)

    poly = np.polyfit(x, pieces, deg=2)
    poly_time = np.polyfit(x, time, deg=2)

    fig, ax = plt.subplots()
    # ax.plot(pieces, label='polynomial pieces')
    ax.plot(x, pieces, label="observation", linewidth=2)
    ax.plot(x, np.polyval(poly, x), label='fit', linewidth=2)
    ax.legend(prop={'size': 16})
    ax.set_ylabel("piece-wise polynomial size", fontdict={'fontsize':16})
    ax.set_xlabel("curve size", fontdict={'fontsize':16})
    s = poly
    # plt.text(80, 4, s, bbox=dict(fill=False, edgecolor='red', linewidth=2))
    plt.savefig("polynomial_pieces_" + dataset)
    
    fig, ax = plt.subplots()
    # ax.plot(time, label='execution time')
    ax.plot(x, time, label="observation")
    ax.plot(x, np.polyval(poly_time, x), label='fit')
    ax.legend(prop={'size': 16})
    ax.set_ylabel("execution time (milliseconds)", fontdict={'fontsize':16})
    ax.set_xlabel("curve size", fontdict={'fontsize':16})
    plt.savefig("time_" + dataset)

plot_results('characters')
plot_results('pigeons')