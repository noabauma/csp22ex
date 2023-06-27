import numpy as np
import matplotlib.pyplot as plt
import platform
import sys

def main():
    

    data = np.loadtxt(sys.argv[1], delimiter=",")
    
    fig, axs = plt.subplots(5)
    fig.suptitle('task1')

    ylabels = [r'$U_L$', r'$\langle  m \rangle$', r'$\chi$', r'$\langle E \rangle$', r'$C_v$']

    for L in np.unique(data[:,0]):
        for i in range(5):
            axs[i].plot(data[data[:,0]==L,1], data[data[:,0]==L,i+2], "x", label="L = " + str(int(L)))
            axs[i].set_xlabel(r'$T$')
            axs[i].set_ylabel(ylabels[i])
            axs[i].legend()

    #print("Critical Temp T_c = ", T_arr[np.argmax(chi_arr)])


    if "WSL" in platform.uname().release or "Microsoft" in platform.uname().release:
        plt.savefig("plots.png")
    else:
        plt.show()

if __name__ == "__main__":
    main()