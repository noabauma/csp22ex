import numpy as np
import matplotlib.pyplot as plt
import platform
import sys

def main():
    data = np.loadtxt(sys.argv[1], delimiter=",")

    L_arr_unique = np.unique(data[:,0]).astype(int)
    n_L = len(L_arr_unique)

    fig, axs = plt.subplots(2,2)
    fig.suptitle('plots')

    for i in range(1,n_L):
        Esys_arr= data[data[:,0]==L_arr_unique[i],1]
        T_arr   = data[data[:,0]==L_arr_unique[i],2]
        M_arr   = data[data[:,0]==L_arr_unique[i],3]
        M_err   = data[data[:,0]==L_arr_unique[i],4]
        chi_arr = data[data[:,0]==L_arr_unique[i],5]
        E_arr   = data[data[:,0]==L_arr_unique[i],6]
        E_err   = data[data[:,0]==L_arr_unique[i],7]
        cv_arr  = data[data[:,0]==L_arr_unique[i],8]


        axs[0,0].errorbar(T_arr, M_arr, yerr=M_err, marker='x', ls='', label="L = " + str(L_arr_unique[i]))

        axs[0,1].plot(T_arr, chi_arr, 'x', label="L = " + str(L_arr_unique[i]))

        axs[1,0].errorbar(T_arr, E_arr, yerr=E_err, marker='x', ls='', label="L = " + str(L_arr_unique[i]))

        axs[1,1].plot(T_arr, cv_arr, 'x', label="L = " + str(L_arr_unique[i]))

        print("Critical Temp T_c = ", T_arr[np.argmax(chi_arr)], " for L = ", L_arr_unique[i])
    
    axs[0,0].set_xlabel(r'$T$')
    axs[0,0].set_ylabel(r'$\langle |M| \rangle$')
    axs[0,0].set_title('magnetization')
    axs[0,0].legend()

    axs[0,1].set_xlabel(r'$T$')
    axs[0,1].set_ylabel(r'$\chi$')
    axs[0,1].set_title('magnetic susceptibility')
    axs[0,1].legend()

    axs[1,0].set_xlabel(r'$T$')
    axs[1,0].set_ylabel(r'$\langle E \rangle$')
    axs[1,0].set_title('energy')
    axs[1,0].legend()

    axs[1,1].set_xlabel(r'$T$')
    axs[1,1].set_ylabel(r'$C_v$')
    axs[1,1].set_title('heat capacity')
    axs[1,1].legend()


    if "WSL" in platform.uname().release or "Microsoft" in platform.uname().release:
        plt.savefig("plots.png")
    else:
        plt.show()


if __name__ == "__main__":
    main()