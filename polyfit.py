import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

truncate = pd.read_csv('trunc_results.csv')
sample = pd.read_csv('sample_results.csv')

pieces_trunc = truncate["number_of_pieces"].to_list()
time_trunc = truncate["time"].to_list()

x = truncate["i"].to_list()

plt.plot(x, pieces_trunc)

poly = np.polyfit(x, pieces_trunc, deg=2)
poly_time = np.polyfit(x, time_trunc, deg=2)

fig, ax = plt.subplots()
# ax.plot(pieces_trunc, label='polynomial pieces')
ax.plot(x, pieces_trunc, label="observation", linewidth=2)
ax.plot(x, np.polyval(poly, x), label='fit', linewidth=2)
ax.legend(prop={'size': 16})
ax.set_ylabel("piece-wise polynomial size", fontdict={'fontsize':16})
ax.set_xlabel("curve size", fontdict={'fontsize':16})
s = poly
# plt.text(80, 4, s, bbox=dict(fill=False, edgecolor='red', linewidth=2))
plt.savefig("polynomial_pieces_truncate_pigeon")
print(poly)
fig, ax = plt.subplots()
# ax.plot(pieces_trunc, label='polynomial pieces')
ax.plot(x, time_trunc, label="observation")
ax.plot(x, np.polyval(poly_time, x), label='fit')
ax.legend(prop={'size': 16})
ax.set_ylabel("execution time (milliseconds)", fontdict={'fontsize':16})
ax.set_xlabel("curve size", fontdict={'fontsize':16})
plt.savefig("time_truncate_pigeon")

print(poly_time)