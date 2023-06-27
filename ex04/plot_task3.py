import numpy as np
import matplotlib.pyplot as plt
import platform
import glob
import pandas as pd

def main():
    filenames = glob.glob("outputs/*.csv")
    filenames = np.sort(filenames)
    n = len(filenames)

    df_info = pd.DataFrame({'algorithm': [], 'L' : [], 'T': [], 'runtime': [], 'tau': [],'MC_speed': []})
    
    L = 32 #we want only this lattice size
    fig, axs = plt.subplots(2,3)
    fig.suptitle("Linear correlation for L = " + str(L))

    num_plot = 0
    for file in filenames:
        info = np.loadtxt(file, dtype=float, delimiter=',', max_rows=1) #L, T, runtime
        data = np.loadtxt(file, dtype=float, delimiter=',', skiprows=1) #nonlinear, linear

        if info[0] == 2:    #L == 2 is trash
            continue

        t = data.shape[1]

        upper_bound = 0
        while upper_bound < t and data[1,upper_bound] > 0.05:
            upper_bound += 1

        m,b = np.polyfit(np.arange(upper_bound),np.log(data[1,:upper_bound]),1)
        tau = -1./m

        title = "M(RT)^2" if "mr2t2" in file else "Wolff"

        arr = np.hstack((title,info,tau,1.0/(info[2]*tau)))
        df_info = df_info.append(pd.DataFrame(arr.reshape(1,-1), columns=list(df_info)), ignore_index=True)
        #df_info = df_info.append(pd.DataFrame(arr).T, ignore_index=True)

        #only plot the ones that are of size L
        if info[0] == L:
            row = int(num_plot/3)
            col = num_plot%3

            num_plot += 1
            T = np.arange(upper_bound)
            axs[row,col].plot(T, data[1,:upper_bound], label=r'$\Phi_{\sigma}(t)$')

            
            axs[row,col].plot(np.exp(m*np.arange(upper_bound)+b), label=fr'fit at $\tau={tau:.0f}$', ls='dashed')

            if col == 0:
                axs[row,col].set_ylabel('Correlation to initial state')
            if row == 1:
                axs[row,col].set_xlabel(r'$t$')

            title = r'$M(RT)^2$' if "mr2t2" in file else r'$Wolff$'

            title += " with Temp = " + str(info[1])

            title += r'$,\tau \approx$' + str(round(tau/L**3,2))  + " sweeps"
            
            axs[row,col].set_title(title)

            axs[row,col].set_ylim((0.0, 1.1))
            axs[row,col].legend()

    print(df_info)

    df_info.to_csv("task3_info.csv")

    if ("Microsoft" in platform.uname().release) or ("WSL" in platform.uname().release):
        plt.savefig("tmp_data/plots.png")
    else:
        plt.show()

if __name__ == "__main__":
    main()