import matplotlib
import pandas 


def save_bar_plot_to_file(directory,log):
        params = directory.split('_')
        l = int(params[2])
        k = int(params[3])
        #print(directory)
        data = pandas.read_csv(directory + '/k_med_scores.csv')
        data["DBA"] = data["DBA"] / data["CDBA"]
        data["WEDGE"] = data["WEDGE"] / data["CDBA"]
        data["CDBA"] = data["CDBA"] / data["CDBA"]
        data = data.drop(columns=["FSA"])
        print(data)
        data = data.set_index("characters")
        ax = data.plot.bar(log=log)
        ax=data.plot(kind='bar')
        ax.legend(loc="best")
        #ax.set_ylim(bottom=data.min().min())
        #print(filename)
        ax.set_ylabel("(k, l)-medians cost")
        ax.get_figure().savefig(directory + '/figures/k_med_scores.png', bbox_inches='tight')
        return 0

