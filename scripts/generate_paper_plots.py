import os
from plot_clustering import plot_clusterings
from initial_clustering_plots import export_initial_clustering_plots
from make_bar_plot import save_bar_plot_to_file

char_exp_dirs = []
pigeon_exp_dir = "results/pigeon_visual_experiment_v_2"
init_exp_dir = "results/initial_clustering_experiment"
movebank_exp_dir = "results/movebank"

#export_initial_clustering_plots(init_exp_dir)
plot_clusterings(pigeon_exp_dir, 25, 600)
plot_clusterings(movebank_exp_dir, 25, 600)
