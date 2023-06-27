#!/usr/bin/python3

from matplotlib import pyplot as plt
import numpy as np
import json
import glob
import re
from sys import argv
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

save = True

datafolder = argv[1] if len(argv) > 1 else "outputs"

folder = ""
if len(argv) < 2:
    folder = "outputs"
else:
    folder = argv[1]

# Find the length of the simulated clusters.
#file = json.load(open("input.json"))
#Ls = file["Ls"]
files = glob.glob(folder + "/magnetization_L*" + ".txt")
Ls = []
for file in files:
    L = re.search("[0-9]+", file).group(0)
    Ls.append(L)


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


def sort(array, col):
    return array[array[:, col].argsort()]


num_Ls = len(Ls)
Ls = np.sort([int(Ls[i]) for i in range(len(Ls))])
color_id = np.linspace(0.3, 1, num_Ls)
fig, axs = plt.subplots(2, 2, figsize=(9.7, 6), sharex=True)

for i, L in enumerate(Ls):
    cid = color_id[i]

    mag_data = np.loadtxt(datafolder + "/magnetization_L"+str(L)+".txt")
    en_data = np.loadtxt(datafolder + "/energies_L"+str(L)+".txt")
    betas = mag_data[:, 0]
    sus = []
    mag = []
    en = []
    cv = []

    for id, beta in enumerate(betas):
	# Note: the c code outputs a magnetization, here we need the total, hence we
	#       divide by the volume.
        M = np.abs(mag_data[id, 1:])
        E = en_data[id, 1:]

        def energy(e):
            return e.mean()

        def magnetization(m):
            return m.mean()

        def specific_heat(e):
            e2 = (e**2).mean()
            e = e.mean()
            return beta**2*(e2-e**2) / L**3

        def susceptibility(m):
            m2 = (m**2).mean()
            m = m.mean()
            return beta*(m2-m**2) / L**3

        sus_val, sus_err = jackKnife(susceptibility, M, 20)
        sus.append([beta, sus_val, sus_err])

        cv_val, cv_err = jackKnife(specific_heat, E, 20)
        cv.append([beta, cv_val, cv_err])

        mag_val, mag_err = jackKnife(magnetization, M, 20)
        mag.append([beta, mag_val, mag_err])

        en_val, en_err = jackKnife(energy, E, 20)
        en.append([beta, en_val, en_err])

    sus = sort(np.array(sus), 0)
    cv = sort(np.array(cv), 0)
    mag = sort(np.array(mag), 0)
    en = sort(np.array(en), 0)

    axs[0, 0].errorbar(mag[:, 0], mag[:, 1]/(L**3), fmt=".-", yerr=mag[:, 2]/(L**3),
                       label=r"$L=$"+str(L), color=plt.cm.PuBu(cid))

    axs[1,0].errorbar(sus[:, 0], sus[:, 1], fmt=".-", yerr=sus[:, 2],
                label=r"$L=$"+str(L), color=plt.cm.PuBu(cid))

    if i == 0:
        axins = inset_axes(axs[0,1], width=1.5, height=1.2, loc='lower right',
                           bbox_to_anchor=(1-0.45, 1-0.7, .3, .3), bbox_transform=axs[1, 0].transAxes)
    axins.semilogy(sus[:, 0], sus[:, 1], '.-', markersize=2,
                   linewidth=1.3, color=plt.cm.PuBu(cid))
    axins.set_xlim(0.65, 0.72)
    axins.set_ylim(0.1, np.max(sus[:, 1])*1.05)
    axins.yaxis.tick_right()
    mark_inset(axs[1, 0], axins, loc1=1, loc2=4, ec='gray', alpha=0.1)

    axs[0, 1].errorbar(en[:, 0], en[:, 1]/(L**3), fmt=".-", yerr=en[:, 2]/(L**3),
                       label=r"$L=$"+str(L), color=plt.cm.PuBu(cid))
    
    axs[1, 1].errorbar(cv[:, 0], cv[:, 1], fmt=".-", yerr=cv[:, 2],
                       label=r"$L=$"+str(L), color=plt.cm.PuBu(cid))


axs[0, 0].set_ylabel(r"$|\vec{m}|$")
axs[0, 0].text(0.05, 0.95, "a", transform=axs[0, 0].transAxes,
               fontweight='bold', va='top')


axs[1, 0].set_ylabel(r"$\chi$")
axs[1, 0].text(0.05, 0.95, "b", transform=axs[1, 0].transAxes,
               fontweight='bold', va='top')

axs[0, 1].set_ylabel(r"$E$")
axs[0, 1].text(0.05, 0.05, "c", transform=axs[0, 1].transAxes,
               fontweight='bold', va='bottom')

axs[1, 1].set_xlabel(r"$1/T$")
axs[1, 1].set_ylabel(r"$c_V$")
axs[1, 1].legend(loc="best")
axs[1, 1].text(0.05, 0.95, "d", transform=axs[1, 1].transAxes,
               fontweight='bold', va='top')

if save:
    plt.savefig('outputs/thdyn.pdf', format='pdf', bbox_inches = "tight")
if not save:
    plt.show()
