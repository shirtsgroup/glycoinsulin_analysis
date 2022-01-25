import glob
import pickle
import natsort
import scipy.stats
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
from matplotlib.path import Path
from MDAnalysis.analysis import dihedrals
from MDAnalysis.analysis.data.filenames import Rama_ref

def read_experimental_data():
    f = open('../Exp_data/experimental_data.txt')
    lines = f.readlines()
    f.close()

    data, height, err = [], [], []
    for i in lines:
        data.append(float(i.split(',')[-1]))

    for i in range(len(data)):
        if i % 2 == 0:   # height
            height.append(data[i])
        else:
            err.append(data[i] - data[i - 1])

    return height, err

if __name__ == "__main__":
    rc("font", **{"family": "sans-serif", "sans-serif": ["DejaVu Sans"], "size": 10})
    # Set the font used for MathJax - more on this later
    rc("mathtext", **{"default": "regular"})
    plt.rc("font", family="Arial")

    folder = ["4EYD", "4EY9", "4EY1", "3I3Z", "2MVC"]
    sys = ["WT"]
    sys.extend([f"GF {i}" for i in range(2, 14)])

    ncolors = 3  # number of distinct colors
    first_rgb = 80  # first value in rgb range
    last_rgb = 10
    cmap = plt.cm.spring  # or whatever color map you choose from matplotlib
    #colors = np.array([cmap(i) for i in np.linspace(first_rgb, last_rgb, ncolors).astype(int)])
    colors = np.array([cmap(i) for i in [40, 110, 180]])

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

    #upper_bounds = np.array([0.025, 100, 37, 101.5])    # upper bounds in x
    #lower_bounds = np.array([-0.006, -12, -5, 86])    # lower bounds in x
    upper_bounds = np.array([115, 115, 115, 115])
    lower_bounds = np.array([-15, -15, -15, -15])
    ranges = upper_bounds - lower_bounds
    annotate_x = [74, 74, 74, 74]  # [0.0163, 67, 24.5, 97.1]
    adjust_x = [3, 3, 3, 3]  # [0.001, 3, 1, 0.25]
    label_x = [-0.0055, -10, -4.5, 86.2]

    # Plot the correlation plots
    exp_h, exp_err = read_experimental_data()
    avg, err = np.mean(beta_all, axis=0), np.std(beta_all, axis=0)
    c_list, e_list = [], []  # Pearson correlation coefficients and the p-values
    for i in range(4):
        c, e = scipy.stats.pearsonr(avg[i], exp_h)
        c_list.append(c)
        e_list.append(e)
    
    plt.figure(figsize=(10, 8))
    for i in range(4):    # 4 examined sites
        plt.subplot(2, 2, i + 1)
        with open(f'beta_data_B{i + 22}.pickle', 'wb') as handle:
            pickle.dump([avg[i], err[i]], handle, protocol=pickle.HIGHEST_PROTOCOL)
        for j in range(len(avg[i])):
            if j == 0:
                plt.errorbar(avg[i][j], exp_h[j], xerr=err[i][j], yerr=exp_err[j], fmt="o", color='red', capsize=2, label='WT')
                plt.text(avg[i][j] + adjust_x[i], exp_h[j] + 0.3, sys[j], weight='bold')
                
            else:
                if j == 1:
                    plt.errorbar(avg[i][j], exp_h[j], xerr=err[i][j], yerr=exp_err[j], fmt="o", color='blue', capsize=2, label='GF')
                else:
                    plt.errorbar(avg[i][j], exp_h[j], xerr=err[i][j], yerr=exp_err[j], fmt="o", color='blue', capsize=2)
                pass_list_1 = [2, 3, 4, 5, 6]
                pass_list_2 = [3, 6]
                pass_list_3 = [3, 4, 6]
                pass_list_4 = [3, 5, 6]
                if i == 0 and j in pass_list_1:
                    pass
                elif i == 0 and (j == 7):
                    plt.text(avg[i][j] + adjust_x[i], exp_h[j] + 0.2, 'GF 3, GF 6, GF 8')
                elif i == 0 and (j == 10):
                    plt.text(avg[i][j] + adjust_x[i], exp_h[j] + 0.2, 'GF 4, GF 7, GF 11')
                elif i == 0 and (j == 8):
                    plt.text(avg[i][j] + adjust_x[i], exp_h[j] + 0.5, 'GF 9')
                elif i == 0 and (j == 11):
                    plt.text(avg[i][j] + adjust_x[i], exp_h[j] - 0.5, 'GF 5, GF 12')


                elif i == 1 and j in pass_list_2:
                    pass
                elif i == 1 and (j == 2):
                    plt.text(avg[i][j] + adjust_x[i], exp_h[j] - 0.8, 'GF 3')
                elif i == 1 and (j == 5):
                    plt.text(avg[i][j] + adjust_x[i], exp_h[j] + 0.5, 'GF 6')
                elif i == 1 and (j == 10):
                    plt.text(avg[i][j] + adjust_x[i], exp_h[j] + 0.2, 'GF 4, GF 7, GF 11')
                
                elif i == 2 and j in pass_list_3:
                    pass
                elif i == 2 and (j == 2):
                    plt.text(avg[i][j] + adjust_x[i], exp_h[j] - 0.8, 'GF 3')
                elif i == 2 and (j == 8):
                    plt.text(avg[i][j] + adjust_x[i], exp_h[j] + 0.5, 'GF 9')
                elif i == 2 and (j == 10):
                    plt.text(avg[i][j] + adjust_x[i], exp_h[j] + 0.2, 'GF 4, GF 7, GF 11')
                elif i == 2 and (j == 11):
                    plt.text(avg[i][j] + adjust_x[i], exp_h[j] + 0.5, 'GF 5, GF 12')

                elif i == 3 and j in pass_list_4:
                    pass
                elif i == 3 and (j == 4):
                    plt.text(avg[i][j] - 4.5 * adjust_x[i], exp_h[j] + 0.2, 'GF 5')
                elif i == 3 and (j == 8):
                    plt.text(avg[i][j] + adjust_x[i], exp_h[j] + 0.5, 'GF 9')
                elif i == 3 and (j == 11):
                    plt.text(avg[i][j] + adjust_x[i], exp_h[j] - 0.8, 'GF 12')
                elif i == 3 and (j == 10):
                    plt.text(avg[i][j] - 10 * adjust_x[i], exp_h[j] + 0.2, 'GF 4, GF 7, GF 11')
                elif i == 3 and (j == 2):
                    plt.text(avg[i][j] + adjust_x[i], exp_h[j] - 0.8, 'GF 3')
                elif i == 3 and (j == 7):
                    plt.text(avg[i][j] - 2 * adjust_x[i], exp_h[j] + 0.2, 'GF 6, GF 8')

                else:
                    plt.text(avg[i][j] + adjust_x[i], exp_h[j] + 0.2, sys[j])

        plt.xlabel(r'Avg. $\beta$-sheet propensity of residue ' + f'B{i + 22} (%)', size=12)
        plt.ylabel(r'$\alpha$-chymotrypsin half-life (min)', size=12)
        plt.grid()    
    
        plt.fill_between([lower_bounds[i], upper_bounds[i]], 2, 8, color=colors[0], alpha=0.3, zorder=0)
        plt.fill_between([lower_bounds[i], upper_bounds[i]], 8, 12.5, color=colors[1], alpha=0.3, zorder=0)
        plt.fill_between([lower_bounds[i], upper_bounds[i]], 12.5, 26, color=colors[2], alpha=0.3, zorder=0)
        #plt.text(label_x[i], 12.8, '(Longer)')
        #plt.text(label_x[i], 8.3, '(Comparable)')
        #plt.text(label_x[i], 2.8, '(Shorter)')
        plt.text(annotate_x[i], 24.6, r'($\tau=$' + f'{c_list[i]:.3f}, p={e_list[i]:.3f})')
        
        plt.xlim([lower_bounds[i], upper_bounds[i]])
        plt.ylim([2, 26])

        if i == 3:
            plt.legend(bbox_to_anchor=(0.23, 0.92))
        else:
            plt.legend(bbox_to_anchor=(0.98, 0.92))
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    #plt.suptitle(r'Correlation b/w $\alpha$-chymotrypsin half-life and $\beta$-sheet propensity of P1$-$P3 region', weight='bold', size=15)
    plt.savefig('avg_beta_propensity_correlation.png', dpi=600)
    
    
    