import time
import pickle
import scipy.stats
import itertools
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib import rc 

def read_pickled_data(path):
    with open(path, 'rb') as handle:
        data = pickle.load(handle)
    return data

def correlation_bootstrapping(x, y, x_err, y_err, n_boot=500):
    r_list = []
    for i in range(n_boot):
        data_1, data_2 = [], []
        for j in range(len(x)):
            data_1.append(np.random.normal(x[j], x_err[j]))
            data_2.append(np.random.normal(y[j], y_err[j]))
        coef, p_val = scipy.stats.pearsonr(data_1, data_2)
        r_list.append(coef)
    
    r = np.mean(r_list)
    r_err = np.std(r_list)

    return r, r_err

if __name__ == "__main__":
    t1 = time.time()
    # Note that this code can only works when the pickled data are present in the folders Metric_1, Metric_2, and Metric_3
    
    # Step 0: Setting things up
    rc("font", **{"family": "sans-serif", "sans-serif": ["DejaVu Sans"], "size": 10})
    # Set the font used for MathJax - more on this later
    rc("mathtext", **{"default": "regular"})
    plt.rc("font", family="Arial")
    ncolors = 13  # number of distinct colors
    first_rgb = 50  # first value in rgb range
    last_rgb = 225
    cmap = plt.cm.jet  # or whatever color map you choose from matplotlib
    colors = np.array([cmap(i) for i in np.linspace(first_rgb, last_rgb, ncolors).astype(int)])

    labels = ['WT']
    labels.extend([f'GF {i}' for i in range(2, 14)])

    # Step 1: Load in the data
    var_list = []
    var_list.append(['SASA of the scissile bond between B25 and B26', 'sasa_sc_B25'])
    var_list.append(['SASA of the scissile bond between B26 and B27', 'sasa_sc_B26'])
    var_list.append(['SASA of the residue B24', 'sasa_res_B24'])
    var_list.append(['SASA of the residue B25', 'sasa_res_B25'])
    var_list.append([r'$\beta$-sheet propensity of B22', 'beta_B22'])
    var_list.append([r'$\beta$-sheet propensity of B23', 'beta_B23'])
    var_list.append([r'$\beta$-sheet propensity of B24', 'beta_B24'])
    var_list.append([r'$\beta$-sheet propensity of B25', 'beta_B25'])

    # Note that the pickled SASA data for each site is based on 5 WT models, with the shape as (5, 13).
    sasa_pb = read_pickled_data('../../proteolytic_stability/Metric_1/sasa_data.pickle')
    A, A_err, B, B_err = np.mean(sasa_pb[0], axis=0), np.linalg.norm(sasa_pb[1], axis=0), np.mean(sasa_pb[2], axis=0), np.linalg.norm(sasa_pb[3], axis=0)
    
    sasa_p1 = read_pickled_data('../../proteolytic_stability/Metric_2/sasa_data.pickle')
    C, C_err, D, D_err = np.mean(sasa_p1[0], axis=0), np.linalg.norm(sasa_p1[1], axis=0), np.mean(sasa_p1[2], axis=0), np.linalg.norm(sasa_p1[3], axis=0)
    
    beta_data, beta_err = [], []
    for i in range(4): # beta-sheet propensity for 4 sites
        data = read_pickled_data(f'../../proteolytic_stability/Metric_3/beta_data_B{i + 22}.pickle')
        beta_data.append(data[0])
        beta_err.append(data[1])
    
    data_all = [A, B, C, D]
    data_all.extend(beta_data)

    err_all = [A_err, B_err, C_err, D_err]
    err_all.extend(beta_err)
    
    # Step 2: Plot the figures
    pairs = list(itertools.combinations(range(8), 2))
    for p in pairs:
        print(f'\nPlotting {var_list[p[0]][0]} against {var_list[p[1]][0]} ...')
        # coef, p_val = scipy.stats.kendalltau(data_all[p[0]], data_all[p[1]])
        r, r_err = correlation_bootstrapping(data_all[p[0]], data_all[p[1]], err_all[p[0]], err_all[p[1]])
        plt.figure()
        for i in range(13):
            plt.errorbar(data_all[p[0]][i], data_all[p[1]][i], xerr=err_all[p[0]][i], yerr=err_all[p[1]][i], fmt="o", capsize=2, color=colors[i], label=labels[i])
        #plt.text(0.02, 0.95, r'(r=' + f'{coef:.3f}, p={p_val:.3f})', transform=plt.gca().transAxes)
        plt.text(0.02, 0.95, f'r={r:.3f} $\pm$ {r_err:.3f}', transform=plt.gca().transAxes)
        plt.xlabel(var_list[p[0]][0])
        plt.ylabel(var_list[p[1]][0])
        plt.grid()

        plt.tight_layout()
        plt.legend()
        plt.savefig(f'figures/{var_list[p[0]][1]}__{var_list[p[1]][1]}.png', dpi=600)

    t2 = time.time()
    print(f'\nTime elapsed: {t2 - t1:.1f} seconds.')

export CMAKE_LIBRARY_PATH="/home/wei-tse/Documents/Software/PLUMED/plumed2_build/:$CMAKE_LIBRARY_PATH"