import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import platform
import sys

def func(x, tau):
    return np.exp(-x/tau)

def main():
    #load data
    data = np.loadtxt(sys.argv[1], delimiter=",")

    temp = sys.argv[2]

    #print("shape = ", data.shape)
    
    times  = np.arange(data.shape[1])
    Phi    = data[0,:]
    Phi_nl = data[1,:]

    popt1, pcov1 = curve_fit(func, times, Phi)
    popt2, pcov2 = curve_fit(func, times, Phi_nl)

    tau    = popt1[0]
    tau_nl = popt2[0]

    print("tau = ", tau)
    print("tau_nl = ", tau_nl)

    fig, axs = plt.subplots(2)
    fig.suptitle('correlation functions at T = ' + str(temp))

    axs[0].plot(times, Phi, label=r'$\Phi(t)$')
    axs[0].plot(times, np.exp(-times/tau), label=r'$\exp({-\frac{t}{\tau}})$')
    axs[0].legend()
    axs[0].set_xlabel(r'time $t$')
    axs[0].set_ylabel(r'$\Phi(t)$')
    axs[0].set_title(r'linear spin-spin correlation function $\tau = $' + str(tau))

    axs[1].plot(times, Phi_nl, label=r'$\Phi_{nl}(t)$')
    axs[1].plot(times, np.exp(-times/tau_nl), label=r'$\exp({-\frac{t}{\tau_{nl}}})$')
    axs[1].legend()
    axs[1].set_xlabel(r'time $t$')
    axs[1].set_ylabel(r'$\Phi_{nl}(t)$')
    axs[1].set_title(r'non-linear spin-spin correlation function $\tau_{nl} = $' + str(tau_nl))


    if ("Microsoft" in platform.uname().release) or ("WSL" in platform.uname().release):
        plt.savefig("tmp_data/plots.png")
    else:
        plt.show()


if __name__ == "__main__":
    main()