import numpy as np
import matplotlib.pyplot as plt
import platform
import sys

def main():
    #load data
    data = np.loadtxt(sys.argv[1], delimiter=",")

    temps = data[:,0]
    M_arr = data[:,1]
    M_err = data[:,2]
    chi_arr = data[:,3]
    E_arr = data[:,4]
    E_err = data[:,5]
    cv_arr = data[:,6]

    fig, axs = plt.subplots(2,2)
    fig.suptitle('plots')

    axs[0,0].errorbar(temps, M_arr, yerr=M_err)
    axs[0,0].set_xlabel(r'$T$')
    axs[0,0].set_ylabel(r'$\langle |M| \rangle$')
    axs[0,0].set_title('magnetization')

    axs[0,1].plot(temps, chi_arr)
    axs[0,1].set_xlabel(r'$T$')
    axs[0,1].set_ylabel(r'$\chi$')
    axs[0,1].set_title('magnetic susceptibility')

    axs[1,0].errorbar(temps, E_arr, yerr=E_err)
    axs[1,0].set_xlabel(r'$T$')
    axs[1,0].set_ylabel(r'$\langle E \rangle$')
    axs[1,0].set_title('energy')

    axs[1,1].plot(temps, cv_arr)
    axs[1,1].set_xlabel(r'$T$')
    axs[1,1].set_ylabel(r'$C_v$')
    axs[1,1].set_title('heat capacity')


    if "WSL" or "Microsoft" in platform.uname().release:
        plt.savefig("plots.png")
    else:
        plt.show()


if __name__ == "__main__":
    main()