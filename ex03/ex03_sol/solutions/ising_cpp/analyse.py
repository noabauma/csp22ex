from matplotlib import pyplot as plt
import numpy as np
import sys
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

beta_c = np.log(1+np.sqrt(2))/2

def exact_M(beta, J=1):
	if beta < beta_c:
		return 0
	else:
		return (1-np.sinh(2*beta*J)**(-4))**(0.125)

beta_cont = np.linspace(0, 1, 100)

dim = int(sys.argv[1])

num_Ls = len(sys.argv)-2  # number of system sizes simulated

fig, axs = plt.subplots(1, 2, figsize=(10, 4), sharex=True)
color_id = np.linspace(0.3, 1, num_Ls)
for i in range(num_Ls):
    L = int(sys.argv[i+2])
    data = np.loadtxt("outputs/thdyn_%id_L%i.txt" % (dim, L))
    
    cid = color_id[i]

    Ed = data[:, 1]
    ids = np.nonzero(Ed)[0].tolist()
    beta = 0.25 * np.log(1. + 4. / Ed[ids])  # 1 / Ed[ids]
    E = data[ids, 2]  # energy density
    M = data[ids, 3]  #  magnetisation density

    axs[0].plot(beta, M, '.-', label=r"$L=%i$" % L,
                   color=plt.cm.PuBu(cid), linewidth=0.7)
    if dim == 2:
        axs[0].plot(beta_cont, [exact_M(b) for b in beta_cont], color='salmon', alpha=0.5)
    axs[0].set_ylabel(r"$|M|/N$")
    axs[0].text(0.05, 0.95, "a", transform=axs[0].transAxes,
                fontweight='bold', va='top')
    axs[0].set_xlim(0.9*np.min(beta), 1.05*np.max(beta))

    axs[1].plot(beta, E, '.-', label=r"$L=%i$" % L,
                   color=plt.cm.PuBu(cid), linewidth=0.7)
    axs[1].set_ylabel(r"$E/N$ $[J]$")
    axs[1].text(0.05, 0.05, "b", transform=axs[1].transAxes,
                fontweight='bold', va='bottom')
    axs[1].set_xlabel(r"$\beta$ $[J/k_{\rm B}]$")
    axs[1].legend(loc='upper right')

plt.savefig('outputs/creutz_thdyn_ising%id.pdf' %dim, bbox_inches="tight")
plt.show()
