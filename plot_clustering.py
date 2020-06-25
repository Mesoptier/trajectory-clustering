import sys
import numpy as np
import matplotlib.pyplot as plt

def plot_clustering(file, figure_name, radius):
    
    colors = ['m', 'g', 'c', 'r', 'y', 'b']
    ax = plt.gca()
    k = int(file.readline())
    clustering = []
    for i in range(k):
        clustering.append({
            "curves": [],
            "center": []
        })
        number_of_curves = int(file.readline())

        for j in range(number_of_curves):
            clustering[-1]["curves"].append(
                    list(map(float, file.readline().rstrip().split(" ")))
                )

            if len(clustering[-1]["curves"][-1]) % 2 != 0:
                print("error: A point on a curve is missing a coordinate")
                sys.exit(0)

        clustering[-1]["center"] = list(map(float, file.readline().rstrip().split(" ")))



    for i in range(len(clustering)):
        cluster = clustering[i]
        for curve in cluster["curves"]:
            for j in range(0, len(curve)-2, 2):
                plt.plot((curve[j], curve[j+2]), (curve[j+1], curve[j+3]), colors[i], linewidth=0.25)

    for i in range(len(clustering)):
        cluster = clustering[i]
        center = cluster["center"]
        for j in range(0, len(center)-2, 2):
            plt.plot((center[j], center[j+2]), (center[j+1], center[j+3]), colors[i], linewidth=3)
            ax.add_artist(plt.Circle((center[j], center[j+1]), radius, color=colors[i]))
        ax.add_artist(plt.Circle((center[-2], center[-1]), radius, color=colors[i]))
  

    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ax.set_aspect('equal', adjustable='box')
    ax.axis('off')
    plt.savefig('figures/' + str(figure_name) + '.png')


    




if __name__ == "__main__":
    
    filename = sys.argv[1]
    figure_name = sys.argv[2]
    radius = float(sys.argv[3])
    file = open(filename, "r")
    plot_clustering(file, figure_name, radius)
    
    
