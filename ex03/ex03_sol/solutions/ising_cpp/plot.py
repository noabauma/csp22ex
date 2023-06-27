#!/usr/bin/python3

from matplotlib import pyplot as plt
import numpy as np
import sys
from mpl_toolkits.axes_grid1.inset_locator import inset_axes,mark_inset


def exact_M(beta, J=1):
	if beta<beta_c: 
		return 0
	else:
		return (1-np.sinh(2*beta*J)**(-4))**(0.125)

beta_cont = np.linspace(0, 1, 100)

fig, axs = plt.subplots(2,2, figsize=(8,6), sharex=True)

dim = int(sys.argv[1])
if dim == 2:
	beta_c = np.log(1+np.sqrt(2))/2
elif dim == 3:
	beta_c = 1/4.51
else:
	beta_c = None

num_Ls = len(sys.argv)-2 # number of system sizes simulated
color_id = np.linspace(0.3, 1, num_Ls)
for i in range(num_Ls):
	L = int(sys.argv[i+2])
	data = np.loadtxt("outputs/thdyn_%id_L%i.txt" %(dim, L))

	cid = color_id[i]

	beta = data[:,0]
	E = data[:,1] # energy density
	M = data[:,2] # magnetisation density
	E2 = data[:,3]
	M2 = data[:,4]

	c_V = beta**2 * (E2 - E**2) * L**dim # specific heat is the susceptibility of total energy per spin
	chi = beta * (M2 - M**2) * L**dim # susceptibility of total magnetisation per spin

	print("The estimate for the critical temperature is", 1/beta[np.argmax(c_V)])

	# Plot magnetisation density
	axs[0,0].plot(beta, M, '.-', label=r"$L=%i$" %L, 
		color=plt.cm.PuBu(cid), linewidth=0.7)
	if dim == 2:
		axs[0].plot(beta_cont, [exact_M(b) for b in beta_cont], color='salmon', alpha=0.5)
	axs[0,0].set_ylabel(r"$|M|/N$")
	axs[0,0].legend(loc='lower right')
	axs[0,0].text(0.05, 0.95, "a", transform=axs[0,0].transAxes, fontweight='bold', va='top')

	# Plot magnetic susceptibility
	axs[0,1].plot(beta, chi, '.-', label=r"$L=%i$" %L, 
		color=plt.cm.PuBu(cid), linewidth=0.7)
	axs[0,1].set_ylabel(r"$\chi$ $[J/k_{\rm B}]$")
	if dim == 2 or dim ==3:
		if i == 0:
			axins_chi = inset_axes(axs[0,1], width=1, height=0.8, loc='upper right', 
				bbox_to_anchor=(1-0.4,1-0.3,.3,.3), bbox_transform=axs[0,1].transAxes)
		axins_chi.semilogy(beta, chi, color=plt.cm.PuBu(cid))
		axins_chi.axvline(beta_c, color='salmon', linewidth=4, alpha=0.5, zorder=0)
		axins_chi.set_xlim(0.8*beta_c, 1.2*beta_c)
		axins_chi.set_ylim(0.1, np.max(chi)*1.05)
		axins_chi.yaxis.tick_right()
		axins_chi.set_xticks([0.8*beta_c, beta_c, 1.2*beta_c])
		axins_chi.set_xticklabels([round(0.8*beta_c, 2), r"$1/T_c$", round(1.2*beta_c, 2)])
		mark_inset(axs[0,1], axins_chi, loc1=2, loc2=4, ec='gray', alpha=0.1)
	axs[0,1].text(0.05, 0.95, "b", transform=axs[0,1].transAxes, fontweight='bold', va='top')

	# Plot energy density
	axs[1,0].plot(beta, E, '.-', label=r"$L=%i$" %L, 
		color=plt.cm.PuBu(cid), linewidth=0.7)
	axs[1,0].set_ylabel(r"$E/N$ $[J]$")
	axs[1,0].text(0.05, 0.05, "c", transform=axs[1,0].transAxes, fontweight='bold', va='bottom')

	# Plot specific heat
	axs[1,1].plot(beta, c_V, '.-', label=r"$L=%i$" %L, 
		color=plt.cm.PuBu(cid), linewidth=0.7)
	axs[1,1].set_xlabel(r"$1/T$ $[J/k_{\rm B}]$")	
	axs[1,1].set_ylabel(r"$c_V$ $[J^4/k_{\rm B}]$")
	
	if dim == 2 or dim ==3:
		if i == 0:
			axins_cv = inset_axes(axs[1,1], width=1, height=0.8, loc='upper right')
		axins_cv.plot(beta, c_V, color=plt.cm.PuBu(cid))
		axins_cv.axvline(beta_c, color='salmon', linewidth=4, alpha=0.5, zorder=0)
		axins_cv.set_xlim(0.8*beta_c, 1.2*beta_c)
		axins_cv.set_ylim(0.6, np.max(c_V)*1.05)
		axins_cv.yaxis.tick_right()
		axins_cv.set_xticks([0.8*beta_c, beta_c, 1.2*beta_c])
		axins_cv.set_xticklabels([round(0.8*beta_c, 2), r"$1/T_c$", round(1.2*beta_c, 2)])
		mark_inset(axs[1,1], axins_cv, loc1=2, loc2=4, ec='gray', alpha=0.1)
	axs[1,1].text(0.05, 0.95, "d", transform=axs[1,1].transAxes, fontweight='bold', va='top')

fig.tight_layout()
plt.savefig('outputs/thdyn_ising%id.pdf' %dim, bbox_inches = "tight")
plt.show()
