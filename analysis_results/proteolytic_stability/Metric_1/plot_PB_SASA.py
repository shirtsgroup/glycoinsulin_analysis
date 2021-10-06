import os
import time
import glob
import pickle
import natsort
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from pymbar import timeseries


def analyze_sasa(xvg_file):
    data = np.transpose(np.loadtxt(xvg_file, comments=["@", "#"]))
    sasa = data[2]

    print(f"Analyzing {xvg_file} ...\n")
    [t, g, N_eff] = timeseries.detectEquilibration(sasa, nskip=1)
    avg = np.mean(sasa[t:])
    std = np.std(sasa) / np.sqrt(len(sasa) / g)

    return avg, std


def sum_up_data(avg_data, std_data):
    avg_data = np.transpose(np.array(avg_data))
    std_data = np.transpose(np.array(std_data))

    avg_all, std_all = [], []  # avg of all GFs

    for i in range(len(avg_data)):
        avg_all.append(np.mean(avg_data[i]))
        std_all.append(1 / len(avg_data) * np.sqrt(sum(map(lambda x: x * x, std_data[i]))))

    return avg_all, std_all


if __name__ == "__main__":
    t1 = time.time()

    rc("font", **{"family": "sans-serif", "sans-serif": ["DejaVu Sans"], "size": 10})
    # Set the font used for MathJax - more on this later
    rc("mathtext", **{"default": "regular"})
    plt.rc("font", family="Arial")

    folders = ["4EYD", "4EY9", "4EY1", "3I3Z", "2MVC"]
    res = ["B25", "B26"]
    exp = [
        "$=$",
        "$-$",
        "$-$",
        "$-$",
        r"$\approx$",
        "$-$",
        "$-$",
        "$-$",
        r"$\approx$",
        "$+$",
        "$-$",
        r"$\approx$",
        "$+$",
    ]
    sys = ["WT"]
    sys.extend([f"GF {i}" for i in range(2, 14)])
    y_labels = [
        'SASA of B25 relative to WT (nm$^{2}$)',
        'SASA of B26 relative to WT (nm$^{2}$)',
    ]
    labels = [f"{sys[i]}\n{exp[i]}" for i in range(len(sys))]
    colors = ['black']
    for i in range(len(labels)):
        if "$+$" in labels[i]:
            colors.append('blue')
        elif r"$-$" in labels[i]:
            colors.append('red')
        elif r"$\approx$" in labels[i]:
            colors.append('magenta')

    B25_sasa, B25_err, B26_sasa, B26_err = [], [], [], []
    if os.path.isfile(f'sasa_data.pickle') is True:
        print('Found the pickled file of SASA data! Reading the data now ...')
        with open('sasa_data.pickle', 'rb') as handle:
            sasa_data = pickle.load(handle)
            B25_sasa, B25_err, B26_sasa, B26_err = sasa_data[0], sasa_data[1], sasa_data[2], sasa_data[3]
    else:
        for folder in folders:
            for i in range(len(res)):
                files = natsort.natsorted(glob.glob(f"{folder}/*sasa_{res[i]}.xvg"))

                avg_list, std_list = [], []
                for f in files:
                    avg, std = analyze_sasa(f)
                    avg_list.append(avg)
                    std_list.append(std)

                if i == 0:
                    B25_sasa.append(avg_list)
                    B25_err.append(std_list)
                elif i == 1:
                    B26_sasa.append(avg_list)
                    B26_err.append(std_list)
        sasa_data = [B25_sasa, B25_err, B26_sasa, B26_err]
        with open('sasa_data.pickle', 'wb') as handle:
            pickle.dump(sasa_data, handle, protocol=pickle.HIGHEST_PROTOCOL)

    # Plot the figure for each WT basis
    for i in range(len(folders)):
        fig = plt.figure()
        for j in range(len(res)):
            ax = fig.add_subplot(2, 1, j + 1)
            plt.bar(
                [f"{sys[j]}\n{exp[j]}" for j in range(len(sys))],
                sasa_data[j * 2][i],
                yerr=sasa_data[j * 2 + 1][i],
                capsize=5,
                color="lightblue",
            )
            plt.ylabel("Average SASA (nm$^{2}$)")
            plt.xticks(fontsize=7)
            plt.grid()

            for color, tick in zip(colors, ax.xaxis.get_major_ticks()):
                tick.label1.set_color(color)
        
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.suptitle(f"Cleavage sites SASA of {folders[i]}-based glycoforms")
        plt.savefig(f"{folders[i]}_sasa.png", dpi=600)

    
    # Plot the avg SASA (for the paper)
    # Average across all wildtype basis
    B25_avg, B25_std = sum_up_data(B25_sasa, B25_err)
    B26_avg, B26_std = sum_up_data(B26_sasa, B26_err)

    # SASA relative to the WT
    x = np.arange(13)

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(2, 1, 1)
    plt.bar(
        x,
        B25_avg - B25_avg[0],
        width=0.4,
        yerr=np.sqrt(np.power(B25_std, 2) + B25_std[0] ** 2),
        capsize=3,
    )
    plt.ylabel("SASA of B25-B26 scissile bond (nm$^{2}$)", fontsize=14)
    plt.xticks(np.arange(len(labels)), tuple(labels), fontsize=12)
    plt.xlim([x[0] - 0.8, x[-1] + 0.8])
    plt.ylim([-0.09, 0.04])
    plt.fill_between([x[0] - 0.8, x[-1] + 0.8], 0.04, color='cyan', alpha=0.8)
    plt.fill_between([x[0] - 0.8, x[-1] + 0.8], -0.09, color='lightgreen', alpha=0.8)
    plt.grid()
    for color, tick in zip(colors, ax.xaxis.get_major_ticks()):
        tick.label1.set_color(color)

    ax = fig.add_subplot(2, 1, 2)
    plt.bar(
        x,
        B26_avg - B26_avg[0],
        width=0.4,
        yerr=np.sqrt(np.power(B26_std, 2) + B26_std[0] ** 2),
        capsize=3,
    )
    plt.ylabel("SASA of B26-B27 scissile bond (nm$^{2}$)", fontsize=14)
    plt.xticks(np.arange(len(labels)), tuple(labels), fontsize=12)
    plt.xlim([x[0] - 0.8, x[-1] + 0.8])
    plt.ylim([-0.12, 0.04])
    plt.fill_between([x[0] - 0.8, x[-1] + 0.8], 0.04, color='cyan', alpha=0.8)
    plt.fill_between([x[0] - 0.8, x[-1] + 0.8], -0.12, color='lightgreen', alpha=0.8)
    plt.grid()
    for color, tick in zip(colors, ax.xaxis.get_major_ticks()):
        tick.label1.set_color(color)

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.suptitle('Average SASA of scissile bond compared to WT', weight='bold', fontsize=15)
    plt.savefig("sasa_PB_all.png", dpi=600)

    t2 = time.time()
    print(f"Elapsed time: {t2 - t1} seconds.")
    
