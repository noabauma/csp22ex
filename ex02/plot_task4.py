import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import platform
import sys

from numpy import dtype

def main():
    T_c = 4.51

    #load data
    data = np.loadtxt(sys.argv[1], delimiter=",")
    
    L_arr = np.unique(data[:,0]).astype(int)
    T_arr = np.unique(data[:,1])

    #print(data)

    n_L = len(L_arr)
    n_T = len(T_arr)
    n   = data.shape[0]
    
    if False:
        #first plot X(T,L) over L to find gamma/nu
        data = data[data[:, 1].argsort()]   #sort by temperature
        chi_arr   = data[:,2]
        for i in range(n_T):
            T   = T_arr[i]
            if T > 0:
                chi = chi_arr[i*n_L:(i+1)*n_L]
                plt.plot(L_arr, chi, '-x', label='T = ' + str(T))
        
        plt.xlabel(r'$L$')
        plt.ylabel(r'$\chi$')
        plt.legend()

        if ("Microsoft" in platform.uname().release) or ("WSL" in platform.uname().release):
            plt.savefig("tmp_data/plots.png")
        else:
            plt.show()

    else:
        fig, ax = plt.subplots()
        # adjust the main plot to make room for the sliders
        plt.subplots_adjust(left=0.25, bottom=0.25)

        nu_init = 1.0
        gamma_nu_init = 1.0

        data = data[data[:, 0].argsort()]   #sort by L
        chi_arr   = data[:,2]

        lines = [plt.plot()]*n_L

        for i in range(n_L):
            L   = L_arr[i]
            chi = chi_arr[i*n_T:(i+1)*n_T]
            lines[i], = plt.plot(np.abs(T_arr-T_c)*L**(1.0/nu_init), chi*L**(-gamma_nu_init), '-x', label='L = ' + str(L))
        ax.set_xlabel(r'$|T-T_c|L^{1/\nu}$')
        ax.set_ylabel(r'$\chi L^{-\gamma/\nu}$')
        ax.set_title('data collapse')
        plt.legend()


        # Make a horizontal slider to control nu.
        axfreq = plt.axes([0.25, 0.1, 0.65, 0.03])
        nu_slider = Slider(
            ax=axfreq,
            label=r'$\nu$',
            valmin=0.1,
            valmax=10.0,
            valinit=nu_init,
        )

        # Make a vertically oriented slider to control gamma/nu
        axamp = plt.axes([0.1, 0.25, 0.0225, 0.63])
        gamma_nu_slider = Slider(
            ax=axamp,
            label=r'$\gamma/\nu$',
            valmin=0.1,
            valmax=10.0,
            valinit=gamma_nu_init,
            orientation="vertical"
        )
        
        
        # The function to be called anytime a slider's value changes
        def update(val):
            for i in range(n_L):
                L   = L_arr[i]
                chi = chi_arr[i*n_T:(i+1)*n_T]
                lines[i].set_xdata(np.abs(T_arr-T_c)*L**(1.0/nu_slider.val))
                lines[i].set_ydata(chi*L**(-gamma_nu_slider.val))
                fig.canvas.draw_idle()

        # register the update function with each slider
        nu_slider.on_changed(update)
        gamma_nu_slider.on_changed(update)

        if ("Microsoft" in platform.uname().release) or ("WSL" in platform.uname().release):
            plt.savefig("tmp_data/plots.png")
        else:
            plt.show()
    


if __name__ == "__main__":
    main()