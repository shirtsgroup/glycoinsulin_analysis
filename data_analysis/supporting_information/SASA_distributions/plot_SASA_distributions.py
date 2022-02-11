import os
import time
import glob
import pickle
import natsort
import scipy
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import rc
from pymbar import timeseries

if __name__ == "__main__":
    t1 = time.time()

    rc("font", **{"family": "sans-serif", "sans-serif": ["DejaVu Sans"], "size": 10})
    # Set the font used for MathJax - more on this later
    rc("mathtext", **{"default": "regular"})
    plt.rc("font", family="Arial")

    ncolors = 13  # number of distinct colors
    first_rgb = 0  # first value in rgb range
    last_rgb = 225
    cmap = plt.cm.jet  # or whatever color map you choose from matplotlib
    colors = np.array([cmap(i) for i in np.linspace(first_rgb, last_rgb, ncolors).astype(int)])

    # Metric 1: Peptide bond SASA
    folders = ["4EYD", "4EY9", "4EY1", "3I3Z", "2MVC"]
    res = ["B25", "B26"]
    sys = ["WT"]
    sys.extend([f"GF {i}" for i in range(2, 14)])
    x_labels = [
        'SASA of B25-B26 scissile bond (nm$^{2}$)',
        'SASA of B26-B27 scissile bond (nm$^{2}$)',
    ]

    plt.figure(figsize=(6, 8))
    for i in range(len(res)):
        plt.subplot(2, 1, i + 1)
        if i == 0:
            plt.text(-0.05, 1.20, 'A', transform=plt.gca().transAxes, fontsize=28, fontweight='bold', va='top', ha='right')
        for j in range(1, 14):
            sasa_list = []
            for folder in folders:
                if j == 1:
                    f = f"../../proteolytic_stability/Metric_1/{folder}/{folder.lower()}_sasa_{res[i]}.xvg"
                else:
                    f = f"../../proteolytic_stability/Metric_1/{folder}/glycoform_{j}_ACS_sasa_{res[i]}.xvg"
                data = np.transpose(np.loadtxt(f, comments=["@", "#"]))
                sasa_list.extend(data[2])
            if j == 1:
                ax = sns.kdeplot(sasa_list, label='WT', common_norm=False)
                # ax = sns.histplot(sasa_list, bins=100, label='WT', stat='probability', color=colors[j - 1], alpha=1)
            else:
                ax = sns.kdeplot(sasa_list, label=f'GF {j}', common_norm=False)
                # ax = sns.histplot(sasa_list, bins=100, label=f'GF {j}', stat='probability', color=colors[j - 1], alpha=1)
        plt.xlabel(x_labels[i])
        plt.ylabel('Probability')

        plt.legend()
        plt.grid()

    plt.tight_layout()
    plt.savefig('SASA_PB_distribution.png', dpi=600)

    # Metric 2: P1 site SASA
    folders = ["4EYD", "4EY9", "4EY1", "3I3Z", "2MVC"]
    res = ["B24", "B25"]
    sys = ["WT"]
    sys.extend([f"GF {i}" for i in range(2, 14)])
    x_labels = [
        'SASA of residue B24 (nm$^{2}$)',
        'SASA of residue B25 (nm$^{2}$)',
    ]

    plt.figure(figsize=(6, 8))
    for i in range(len(res)):
        plt.subplot(2, 1, i + 1)
        if i == 0:
            plt.text(-0.05, 1.20, 'B', transform=plt.gca().transAxes, fontsize=28, fontweight='bold', va='top', ha='right')
        for j in range(1, 14):
            sasa_list = []
            for folder in folders:
                if j == 1:
                    f = f"../../proteolytic_stability/Metric_2/{folder}/{folder.lower()}_sasa_res_{res[i]}.xvg"
                else:
                    f = f"../../proteolytic_stability/Metric_2/{folder}/glycoform_{j}_ACS_sasa_res_{res[i]}.xvg"
                data = np.transpose(np.loadtxt(f, comments=["@", "#"]))
                sasa_list.extend(data[2])
            if j == 1:
                ax = sns.kdeplot(sasa_list, label='WT', common_norm=False)
            else:
                ax = sns.kdeplot(sasa_list, label=f'GF {j}', common_norm=False)
        plt.xlabel(x_labels[i])
        plt.ylabel('Probability')

        plt.legend()
        plt.grid()

    plt.tight_layout()
    plt.savefig('SASA_res_distribution.png', dpi=600)
    



