import matplotlib.pyplot as plt

def plot(f, c, lbl, marker='-'):
    lines = f.readlines()
    x = [float(line.split()[0]) for line in lines]
    y = [float(line.split()[1]) for line in lines]
    plt.plot(x, y, marker, c=c, lw=.8, label=lbl)

i = 52
fig = plt.figure(dpi=600, frameon=False)
ax = fig.add_subplot(111)
crv_file = "crv" + str(i)
ii_dtw_file = "ii_dtw" + str(i)
best_file = "best" + str(i)
with open(crv_file) as f:
    plot(f, 'k', 'Original', '.')
with open(ii_dtw_file) as f:
    plot(f, '#EE7733', 'Imai--Iri DTW', '-*')
with open(best_file) as f:
    plot(f, '#0077BB', 'Imai--Iri CDTW', '-x')

ax.set_aspect('equal', adjustable='box')
ax.legend(frameon=False,fontsize='x-small')
plt.axis('off')
plt.savefig(crv_file + ".png", dpi=600, bbox_inches='tight')

i = 134
fig = plt.figure(dpi=600, frameon=False)
ax = fig.add_subplot(111)
crv_file = "crv" + str(i)
gr_dtw_file = "gr_dtw" + str(i)
best_file = "best" + str(i)
with open(crv_file) as f:
    plot(f, 'k', 'Original', '.')
with open(gr_dtw_file) as f:
    plot(f, '#EE7733', 'Greedy DTW', '-*')
with open(best_file) as f:
    plot(f, '#0077BB', 'Imai--Iri CDTW', '-x')

ax.set_aspect('equal', adjustable='box')
ax.legend(frameon=False,fontsize='x-small')
plt.axis('off')
plt.savefig(crv_file + ".png", dpi=600, bbox_inches='tight')

i = 154
fig = plt.figure(dpi=600, frameon=False)
ax = fig.add_subplot(111)
crv_file = "crv" + str(i)
gr_cdtw_file = "gr_cdtw" + str(i)
gr_fr_file = "gr_fr" + str(i)
ii_fr_file = "ii_fr" + str(i)
best_file = "best" + str(i)
with open(crv_file) as f:
    plot(f, 'k', 'Original', '.')
# with open(gr_cdtw_file) as f:
#     plot(f, '#EE7733', 'Greedy CDTW', '-o')
with open(gr_fr_file) as f:
    plot(f, '#009988', 'Greedy Frechet', '-|')
# with open(ii_fr_file) as f:
#     plot(f, '#EE3377', 'Imai--Iri Frechet', '-_')
# with open(best_file) as f:
#     plot(f, '#0077BB', 'Imai--Iri CDTW', '-x')

ax.set_aspect('equal', adjustable='box')
ax.legend(frameon=False,fontsize='x-small')
plt.axis('off')
plt.savefig(crv_file + ".png", dpi=600, bbox_inches='tight')
