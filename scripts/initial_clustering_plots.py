import matplotlib.pyplot as plt

def export_initial_clustering_plots(directory):
    results = open(directory + "/results.txt", 'r')
    results = results.read().rstrip().split("\n")
    dataset = [ "Bladon Heath", "Church Hanborough", "Horspath", "Weston" ]
    labels = ["Gonzalez + PAM", "PAM", "Gonzalez"]
    b_h = []
    ch = []
    hors =[]
    west = []
    for line in results:
        line = line.split("\t")
        tgt_list = b_h
        if line[0] == "Bladon & Church route recapping/bladon heath":
            tgt_list = b_h
        if line[0] == "Bladon & Church route recapping/church hanborough":
            tgt_list = ch
        if line[0] == "Horspath":
            tgt_list = hors
        if line[0] == "Weston":
            tgt_list = west

        tgt_list.append((float(line[2]), float(line[3]), float(line[4])))

    tables = [b_h, ch, hors, west]

    for i in range(len(tables)):
        Y_tuples = tables[i]
        export_plot(range(1, len(Y_tuples)+1), Y_tuples, labels, dataset[i], directory + "/figures/" + dataset[i] + ".pdf")

def export_plot(X, Y_tuples, L, title, filename):
    fig,ax = plt.subplots()
    linestyles = ['-', '-.', '--']
    linewidths= [3,2,2]
    for i in range(3):
        Y = [ y[i] for y in Y_tuples ]
        objects = ax.plot(X, Y, linestyle=linestyles[i], alpha=1., label=L[i], linewidth=linewidths[i])

    plt.legend()
    plt.title(title)
    plt.xlabel("k")
    plt.ylabel("cost")
    plt.savefig(filename)
    plt.clf()
