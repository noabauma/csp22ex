import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import platform
import sys

def func(x, gamma):
    T_c = 4.51
    return np.abs(x - T_c)**-gamma

def main():
    #load data
    data = np.loadtxt(sys.argv[1], delimiter=",")

    data = data[data[:, 0].argsort()]

    #print(data)
    
    temps = data[:,0]
    M_arr = data[:,1]
    M_err = data[:,2]
    chi_arr = data[:,3]
    E_arr = data[:,4]
    E_err = data[:,5]
    cv_arr = data[:,6]

    popt, pcov = curve_fit(func, temps, chi_arr)
    print("[gamma] = ", popt)

    
    fig, axs = plt.subplots(2,2)
    fig.suptitle('plots')

    axs[0,0].errorbar(temps, M_arr, yerr=M_err)
    axs[0,0].set_xlabel(r'$T$')
    axs[0,0].set_ylabel(r'$\langle |M| \rangle$')
    axs[0,0].set_title('magnetization')

    axs[0,1].plot(temps, chi_arr)
    axs[0,1].plot(temps, func(temps, popt[0]), label=r'$|T - T_c|^{-\gamma}\ \gamma =$' + str(round(popt[0],3)))
    axs[0,1].set_xlabel(r'$T$')
    axs[0,1].set_ylabel(r'$\chi$')
    axs[0,1].set_title('magnetic susceptibility')
    axs[0,1].legend()

    axs[1,0].errorbar(temps, E_arr, yerr=E_err)
    axs[1,0].set_xlabel(r'$T$')
    axs[1,0].set_ylabel(r'$\langle E \rangle$')
    axs[1,0].set_title('energy')

    axs[1,1].plot(temps, cv_arr)
    axs[1,1].set_xlabel(r'$T$')
    axs[1,1].set_ylabel(r'$C_v$')
    axs[1,1].set_title('heat capacity')

    print("Critical Temp T_c = ", temps[np.argmax(cv_arr)])


    if ("Microsoft" in platform.uname().release) or ("WSL" in platform.uname().release):
        plt.savefig("tmp_data/plots.png")
    else:
        plt.show()

    """
    plt.plot(temps, data[:,7], label="L = " + str(int(16)))
    plt.xlabel(r'$T$')
    plt.ylabel(r'$U_L$')
    plt.legend()
    plt.ylim((0,1))
    plt.show()
    """


if __name__ == "__main__":
    main()