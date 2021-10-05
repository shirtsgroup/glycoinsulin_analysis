import glob
import pickle
import natsort

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
from matplotlib.path import Path
from MDAnalysis.analysis import dihedrals
from MDAnalysis.analysis.data.filenames import Rama_ref

if __name__ == "__main__":
    rc("font", **{"family": "sans-serif", "sans-serif": ["DejaVu Sans"], "size": 10})
    # Set the font used for MathJax - more on this later
    rc("mathtext", **{"default": "regular"})
    plt.rc("font", family="Arial")

    folder = ["4EYD", "4EY9", "4EY1", "3I3Z", "2MVC"]
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
    labels = [f"{sys[i]}\n{exp[i]}" for i in range(len(sys))]
    colors = ['black']
    for i in range(len(labels)):
        if "$+$" in labels[i]:
            colors.append('blue')
        elif r"$-$" in labels[i]:
            colors.append('red')
        elif r"\approx$" in labels[i]:
            colors.append('magenta')

    beta_all = []
    for i in range(len(folder)):
        beta_GF = []
        pickled_files = natsort.natsorted(glob.glob(f"{folder[i]}_Ramachandran_analysis/{folder[i].lower()}_beta_res*.pickle"))
        for j in range(len(pickled_files)):  # different residues
            print(f'Reading {pickled_files[j]} ... \n')
            with open(pickled_files[j], "rb") as handle:
                beta = pickle.load(handle)  # a list of 13 beta fractions
                beta_GF.append(beta)
        beta_all.append(beta_GF)

    fig = plt.figure(figsize=(16, 8))
    avg, err = np.mean(beta_all, axis=0), np.std(beta_all, axis=0)
    for j in range(4):   # 4 different residues
        ax = fig.add_subplot(2, 2, j + 1)

        plt.bar(
            [f"{sys[i]}\n{exp[i]}" for i in range(len(sys))],
            avg[j],
            yerr=err[j], 
            capsize=3,
            color="lightblue",
        )
        for color, tick in zip(colors, ax.xaxis.get_major_ticks()):
            tick.label1.set_color(color)
        plt.xticks(fontsize=12)
        plt.ylabel(r"Avg. $\beta$-sheet propensity of " + f"B{j + 22} (%)", fontsize=15)
        plt.grid()
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig("avg_beta_propensity.png", dpi=600)
