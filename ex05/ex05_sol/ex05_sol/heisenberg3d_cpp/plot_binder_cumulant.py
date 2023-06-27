#!/usr/bin/python3
# Usage: python binder_cumulant.py <data folder. default = outputs>

from matplotlib import pyplot as plt
from sys import argv
import numpy as np
import glob
import re
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset


save = True

folder = ""
if len(argv) < 2: folder = "outputs"
else: folder = argv[1]


def jackKnife(procedure, x, n):
    block_size = int(np.ceil(len(x) / n))
    f = []
    for i in range(n):
        start = int(i*block_size)
        end = int(min((i+1)*block_size, len(x)-1))
        x_loo = np.concatenate([x[:start], x[end:]])
        f.append(procedure(x_loo))

    f_bar = procedure(x)
    f = np.array(f)

    delta = np.sqrt((n-1)/n * np.sum((f-f_bar)**2))
    estimate = n * f_bar - (n-1) * np.mean(f)
    return estimate, delta


def binderCumulant(x2):
    x2_mean = x2.mean()
    x4_mean = (x2**2).mean()
    return 1 - x4_mean/(3*x2_mean**2)


def sort(array, col):
    return array[array[:, col].argsort()]


# Find the length of the simulated clusters.
files = glob.glob(folder + "/magnetization_L*" + ".txt")
Ls = []
for file in files:
    L = re.search("[0-9]+", file).group(0)
    Ls.append(L)

num_Ls = len(Ls)
Ls = np.sort([int(Ls[i]) for i in range(len(Ls))])
color_id = np.linspace(0.3, 1, num_Ls)
fig, ax = plt.subplots(1, figsize=(7, 4))

cumulants = []
for i, L in enumerate(Ls):
    cid = color_id[i]

    data = np.loadtxt(folder + "/magnetization_L"+str(L)+".txt")
    betas = data[:, 0]
    out = []

    for id, beta in enumerate(betas):
        M = data[id, 1:]
        b, db = jackKnife(binderCumulant, M**2, 20)

        out.append([beta, b, db])

    out = sort(np.array(out), 0)
    cumulants.append(out[:, 1])

    np.savetxt(folder + "/binder_cumulants_L" + str(L) + ".txt", out, header='T\tU\tdelta U')
    ax.errorbar(out[:, 0], out[:, 1], fmt=".-", yerr=out[:, 2], 
                label=r"$L=$"+str(L), color=plt.cm.PuBu(cid))

    if i == 0:
        axins = inset_axes(ax, width=1.8, height=1.5, loc='lower right',
            bbox_to_anchor=(1-0.6, 1-0.7, .3, .3), bbox_transform=ax.transAxes)
    axins.plot(out[:, 0], out[:, 1], '.-', markersize=4, linewidth=1.3, color=plt.cm.PuBu(cid))
    axins.set_xlim(0.683, 0.706)
    axins.set_ylim(0.56, 0.668)
    axins.yaxis.tick_right()
    mark_inset(ax, axins, loc1=1, loc2=3, ec='gray', alpha=0.1)

#cumulants = np.array(cumulants)
#var_cumulants = np.std(cumulants, axis=0)
#print(var_cumulants)

ax.legend(loc="best")
ax.set_xlabel(r"$1/T$")
ax.set_ylabel(r"$U_4$")
#ax.set_title("Binder cumulant")

if save : plt.savefig("outputs/binder_cumulant.pdf", format = 'pdf')
else: plt.show()
