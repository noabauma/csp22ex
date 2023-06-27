import sys
import numpy as np
import matplotlib.pyplot as plt


out_dir = 'outputs/'

num_Ls = len(sys.argv)-1 # number of system sizes simulated
color_id = np.linspace(0.3, 1, num_Ls)

fig, axs = plt.subplots(1, 2, figsize=(10, 4), sharex=True)
for i in range(num_Ls):
	L = int(sys.argv[i+1])

	T = np.loadtxt(out_dir+'T%i.txt' %L)
	M = np.loadtxt(out_dir+'M%i.txt' %L)
	E = np.loadtxt(out_dir+'E%i.txt' %L)
	M_err = np.loadtxt(out_dir+'M_err%i.txt' %L)
	E_err = np.loadtxt(out_dir+'E_err%i.txt' %L)


	# magnetization plot
	axs[0].errorbar(T, M, yerr=E_err, fmt='.-', label='L=%i' %L)
	axs[0].set_xlabel(r'$T$')
	axs[0].set_ylabel(r'$M/N$')
	
	# energy plot
	axs[1].errorbar(T, E, yerr=E_err, fmt='.-', label='L=%i' %L)
	axs[1].set_ylabel(r'$E/N$')
	
plt.legend()
plt.show()