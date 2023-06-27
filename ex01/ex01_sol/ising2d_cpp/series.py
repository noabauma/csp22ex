from matplotlib import pyplot as plt
import numpy as np
import sys
import json

with open('input.json') as f:
    specs = json.load(f)
Ts = specs["Ts"]

L = 16
data = np.loadtxt("outputs/magnetizations_L%i.txt" % L)

t = 3
dat = data[len(Ts)-1-t]
steps = np.arange(len(dat))
    
T = Ts[t]

fig, axs = plt.subplots(1, 2, figsize=(10, 4))
#fig.suptitle(r"$L=%i$" % L)

skip = 50
axs[0].plot(steps[::skip], dat[::skip], '.', markersize=0.8,
              label=r"$T=%.2f<T_c$"%T+"\n"+"$L=%i$"%L)
axs[0].legend()
axs[0].set_xlabel(r"$t$")
axs[0].set_ylabel(r"$M(t)$")
axs[0].set_xlim(0, len(dat))
axs[0].set_ylim(-1.1, 1.1)

weights = np.ones_like(dat) / len(dat)
axs[1].hist(dat, weights=weights, bins=np.linspace(-1, 1, 50))
axs[1].set_xlabel(r"M")
axs[1].set_ylabel(r"$p(M)$")

fig.savefig("outputs/flips_L%i_T%.2f.pdf" % (L, T))