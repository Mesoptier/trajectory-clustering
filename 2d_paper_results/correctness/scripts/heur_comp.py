import matplotlib.pyplot as pt
import pandas as pd

pigeons_frame = pd.read_csv("2d_paper_results/correctness/pigeon_results.csv")
characters_frame = pd.read_csv("2d_paper_results/correctness/character_results.csv")

ax = pigeons_frame.plot.bar(rot=0)
ax.set_ylabel("cdtw cost")
ax.set_xlabel("pairs of pigeon trajectories")
fig = ax.get_figure()
fig.savefig('2d_paper_results/correctness/figures/pigeons_comparison.png')

ax = characters_frame.plot.bar(rot=0)
ax.set_ylabel("cdtw cost")
ax.set_xlabel("pairs of character trajectories")
fig = ax.get_figure()
fig.savefig('2d_paper_results/correctness/figures/characters_comparison.png')
