import numpy as np
import matplotlib.pyplot as plt
import platform
import sys

def main():
    

    data = np.loadtxt(sys.argv[1], delimiter=",")
    
    
    for L in np.unique(data[:,0]):
        plt.plot(data[data[:,0]==L,1], data[data[:,0]==L,2], "x", label="L = " + str(int(L)))
        plt.xlabel(r'$T$')
        plt.ylabel(r'$U_L$')
        plt.legend()
        #plt.ylim((0,1))

    #print("Critical Temp T_c = ", T_arr[np.argmax(chi_arr)])


    if "WSL" in platform.uname().release or "Microsoft" in platform.uname().release:
        plt.savefig("plots.png")
    else:
        plt.show()

if __name__ == "__main__":
    main()