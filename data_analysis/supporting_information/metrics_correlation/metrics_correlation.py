import time
import pickle
import random
import scipy.stats
import itertools
import warnings
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib import rc 

warnings.filterwarnings("ignore", category=RuntimeWarning)


def read_pickled_data(path):
    with open(path, 'rb') as handle:
        data = pickle.load(handle)
    return data

def bootstrapping_normal(x, y, x_err, y_err, n_boot=500):
    """
    This function generates bootstrap samples for x and y data by drawing 
    samples from normal distributions centered at the means of x and y data
    to calculate the correlation coefficient and its uncertainty.

    Note that this method is not used in the most updated code but we keep
    here as a record.
    """
    r_list = []
    for i in range(n_boot):
        data_1, data_2 = [], []
        for j in range(len(x)):
            data_1.append(np.random.normal(x[j], x_err[j]))
            data_2.append(np.random.normal(y[j], y_err[j]))
        coef, p_val = scipy.stats.kendalltau(data_1, data_2)
        r_list.append(coef)
    
    r = np.mean(r_list)
    r_err = np.std(r_list)

    return r, r_err

def bootstrapping_for_metrics(x_data, y_data):
    """
    This functions bootstraps over the raw data of variables x and y, which are 
    all computational values to calcualte the uncertainty of the correlation coefficient.

    Parameters
    ----------
    x_data (array-like): 
        The raw data of x. Should be in the shape of (n_variants, n_WTmodels)
    y_data (array-like): 
        The raw data of x. Should be in the shape of (n_variants, n_WTmodels)
    """
    # Boostrap over the data of 5 different WT models
    x_data = np.transpose(x_data)   # should be (13, 5), or (n_variants, n_WTmodels)
    y_data = np.transpose(y_data)
    r_list = []
    for i in range(500):  # number of bootstrap
        xx, yy = [], []   # bootstrap samples for variables x and y
        for j in range(len(x_data)):  # 13 variants:
            xx.append(np.mean(random.choices(x_data[j], k=len(x_data[j]))))
            yy.append(np.mean(random.choices(y_data[j], k=len(y_data[j]))))
        coef, p_val = scipy.stats.kendalltau(xx, yy)
        r_list.append(coef)
    # r = np.mean(r_list)
    r_err = np.std(r_list)  # we only use bootstrapping to calculate the uncertainty here

    return r_err

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

    # Step 1: Load in the mean value data set (for plotting)
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
    sasa_pb = read_pickled_data('../../proteolytic_stability/Metric_1/sasa_data.pickle')  # shape: (4, 5, 13)
    A, A_err, B, B_err = np.mean(sasa_pb[0], axis=0), np.linalg.norm(sasa_pb[1], axis=0), np.mean(sasa_pb[2], axis=0), np.linalg.norm(sasa_pb[3], axis=0)
    
    sasa_p1 = read_pickled_data('../../proteolytic_stability/Metric_2/sasa_data.pickle')  # shape: (4, 5, 13)
    C, C_err, D, D_err = np.mean(sasa_p1[0], axis=0), np.linalg.norm(sasa_p1[1], axis=0), np.mean(sasa_p1[2], axis=0), np.linalg.norm(sasa_p1[3], axis=0)
    
    beta_mean, beta_err = [], []
    for i in range(4): # beta-sheet propensity for 4 sites
        data = read_pickled_data(f'../../proteolytic_stability/Metric_3/beta_data_B{i + 22}.pickle')
        beta_mean.append(data[0])
        beta_err.append(data[1])
    
    mean_all = [A, B, C, D]
    mean_all.extend(beta_mean)

    err_all = [A_err, B_err, C_err, D_err]
    err_all.extend(beta_err)

    # Step 2: Perpare datasets for bootstsrapping
    a, b = sasa_pb[0], sasa_pb[2]   # shape: (5, 13)
    c, d = sasa_p1[0], sasa_p1[2]   # shape: (5, 13)
    beta_all = read_pickled_data('../../proteolytic_stability/Metric_3/beta_data_all.pickle')  # shape: (5, 4, 13)
    beta_data = [beta_all[:, i, :] for i in range(4)]
    data_all = [a, b, c, d]
    data_all.extend(beta_data)

    # Step 3: Plot the figures
    pairs = list(itertools.combinations(range(8), 2))
    for p in pairs:
        print(f'\nPlotting {var_list[p[0]][0]} against {var_list[p[1]][0]} ...')
        r, _ = scipy.stats.pearsonr(mean_all[p[0]], mean_all[p[1]])
        r_err = bootstrapping_for_metrics(data_all[p[0]], data_all[p[1]])
        
        plt.figure()
        for i in range(13):
            plt.errorbar(mean_all[p[0]][i], mean_all[p[1]][i], xerr=err_all[p[0]][i], yerr=err_all[p[1]][i], fmt="o", capsize=2, color=colors[i], label=labels[i])
        #plt.text(0.02, 0.95, r'(r=' + f'{coef:.3f}, p={p_val:.3f})', transform=plt.gca().transAxes)
        plt.text(0.02, 0.95, f'r={r:.3f} $\pm$ {r_err:.3f}', transform=plt.gca().transAxes)
        plt.xlabel(var_list[p[0]][0])
        plt.ylabel(var_list[p[1]][0])
        plt.grid()

        plt.tight_layout()
        plt.legend()

        # Adjust the x and y range
        if 'beta' in f'{var_list[p[0]][1]}':   # x variable
            plt.xlim([-10, 110])
        if 'beta' in f'{var_list[p[1]][1]}':   # y variable
            plt.ylim([-10, 110])

        # Adjust the position of the legend box
        group_1 = [
            'beta_B24__beta_B25.png',
            'sasa_res_B24__beta_B22.png',
            'sasa_res_B24__beta_B25.png',
            'sasa_res_B25__beta_B22.png',
            'sasa_res_B25__beta_B25.png',
            'sasa_res_B24__beta_B22.png',
            'sasa_res_B25__beta_B22.png',
            'sasa_sc_B25__beta_B22.png',
            'sasa_sc_B25__beta_B25.png',
            'sasa_sc_B25__beta_B24.png',
            'sasa_sc_B26__beta_B25.png',
            'sasa_sc_B26__beta_B22.png',
            'sasa_sc_B26__sasa_res_B24.png',
            'sasa_sc_B26__sasa_res_B25.png', 
            'sasa_sc_B25__sasa_res_B25.png', 
            'sasa_sc_B25__sasa_res_B24.png', 
            'sasa_sc_B25__sasa_sc_B26.png',        
        ]

        group_2 = [
            'sasa_res_B24__sasa_res_B25.png',
        ]

        png_name = f'{var_list[p[0]][1]}__{var_list[p[1]][1]}.png'
        if png_name in group_1:
            plt.legend(bbox_to_anchor=(0.001, 0.61), loc='center left')
        elif png_name in group_2:
            plt.legend(bbox_to_anchor=(0.82, 0.33), loc='center left')

        if 'beta' in png_name and 'sasa' in png_name:
            plt.savefig(f'figures/mixed_metrics/{png_name}', dpi=600)
        else:
            plt.savefig(f'figures/{png_name}', dpi=600)        

    t2 = time.time()
    print(f'\nTime elapsed: {t2 - t1:.1f} seconds.')
