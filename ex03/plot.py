import numpy as np
import matplotlib.pyplot as plt
import platform
import sys

def main():
    L = 6
    N = L**3

    py_mode = False if len(sys.argv[1]) > 0 else True 

    if py_mode:
        out_dir = 'outputs/'
        T_arr = np.loadtxt(out_dir+'T'+str(L)+'.txt')
        M_arr = np.loadtxt(out_dir+'M'+str(L)+'.txt')
        E_arr = np.loadtxt(out_dir+'E'+str(L)+'.txt')
        M_err = np.loadtxt(out_dir+'M_err'+str(L)+'.txt')
        E_err = np.loadtxt(out_dir+'E_err'+str(L)+'.txt')

        beta_arr = 1.0/T_arr

        chi_arr = N*beta_arr*M_err**2 # magnetic susceptibility per spin
        cv_arr  = beta_arr**2*E_err**2
    else:
        data = np.loadtxt(sys.argv[1], delimiter=",")
        
        L_arr   = data[:,0]
        Esys_arr= data[:,1]
        T_arr   = data[:,2]
        M_arr   = data[:,3]
        M_err   = data[:,4]
        chi_arr = data[:,5]
        E_arr   = data[:,6]
        E_err   = data[:,7]
        cv_arr  = data[:,8]

    if True:
        plt.plot(Esys_arr, T_arr, "x")
        plt.xlabel(r"$E_{sys}$")
        plt.ylabel("T")
        plt.title("total energies of the microcanonical ensemble & their temperatures")
        plt.show()
    

    fig, axs = plt.subplots(2,2)
    fig.suptitle("Critical Temp $T_c$ = " + str(T_arr[np.argmax(chi_arr)]) + " for L = 16")

    axs[0,0].errorbar(T_arr, M_arr, yerr=M_err, marker='x', ls='')
    axs[0,0].set_xlabel(r'$T$')
    axs[0,0].set_ylabel(r'$\langle |M| \rangle$')
    axs[0,0].set_title('magnetization')

    axs[0,1].plot(T_arr, chi_arr, 'x')
    axs[0,1].set_xlabel(r'$T$')
    axs[0,1].set_ylabel(r'$\chi$')
    axs[0,1].set_title('magnetic susceptibility')

    axs[1,0].errorbar(T_arr, E_arr, yerr=E_err, marker='x', ls='')
    axs[1,0].set_xlabel(r'$T$')
    axs[1,0].set_ylabel(r'$\langle E \rangle$')
    axs[1,0].set_title('energy')

    axs[1,1].plot(T_arr, cv_arr, 'x')
    axs[1,1].set_xlabel(r'$T$')
    axs[1,1].set_ylabel(r'$C_v$')
    axs[1,1].set_title('heat capacity')

    print("Critical Temp T_c = ", T_arr[np.argmax(chi_arr)])


    if "WSL" in platform.uname().release or "Microsoft" in platform.uname().release:
        plt.savefig("plots.png")
    else:
        plt.show()

if __name__ == "__main__":
    main()