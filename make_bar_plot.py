import matplotlib
import pandas

def save_barplot_to_file(filename,log):
	params = filename.split('_')
	l = int(params[2])
	k = int(params[3])
	data = pandas.read_csv(filename + '/k_med_scores.csv')
	data = data.set_index("characters")
	#ax = data.plot.bar(log=log)
	ax=data.plot(kind='bar')
	ax.legend(loc="best")
	if log:
		filename += '_log'
	# #ax.set_ylim(bottom=data.min().min())
	print(filename)
	ax.set_ylabel("(k, l)-medians cost")
	ax.get_figure().savefig(filename+'.png', bbox_inches='tight')
	return 0

for filename in ['char_exp_12_2_50', 'char_exp_9_2_50', 'char_exp_6_2_50']:
	for log in [True,False]:
		save_barplot_to_file(filename=filename, log=log)
