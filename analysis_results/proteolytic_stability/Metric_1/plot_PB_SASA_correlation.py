import os
import time
import glob
import pickle
import natsort
import scipy
import random
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

def bootstrapping_sample(x_data, y_mean, y_err):
    """
    This functions bootstraps over the raw data of x (e.g. computational values)
    and a normal distribution centered at the mean of y (e.g. experimental values)
    to calcualte the uncertainty of the correlation coefficient.

    Parameters
    ----------
    x (array-like): 
        The raw data of x. Should be in the shape of (n_variants, n_WTmodels)
    y_mean (array-like): 
        The experimental data. The length should be 13. 
    y_err (array-like): 
        The uncertainty of the experimental data points. The length should be 13
    """
    # Boostrap over the data of 5 different WT models
    x_data = np.transpose(x_data)   # should be (13, 5), or (n_variants, n_WTmodels)
    r_list = []
    for i in range(500):  # number of bootstrap
        xx, yy = [], []   # bootstrap samples for variables x and y
        for j in range(len(x_data)):  # 13 variants:
            xx.append(np.mean(random.choices(x_data[j], k=len(x_data[j]))))
            yy.append(np.mean(np.random.normal(y_mean[j], y_err[j], len(x_data[j]))))
        coef, p_val = scipy.stats.kendalltau(xx, yy)
        r_list.append(coef)
    # r = np.mean(r_list)
    r_err = np.std(r_list)  # we only use bootstrapping to calculate the uncertainty here

    return r_err

if __name__ == "__main__":
    t1 = time.time()

    rc("font", **{"family": "sans-serif", "sans-serif": ["DejaVu Sans"], "size": 10})
    # Set the font used for MathJax - more on this later
    rc("mathtext", **{"default": "regular"})
    plt.rc("font", family="Arial")

    folders = ["4EYD", "4EY9", "4EY1", "3I3Z", "2MVC"]
    res = ["B25", "B26"]
    sys = ["WT"]
    sys.extend([f"GF {i}" for i in range(2, 14)])

    ncolors = 3  # number of distinct colors
    first_rgb = 80  # first value in rgb range
    last_rgb = 10
    cmap = plt.cm.spring  # or whatever color map you choose from matplotlib
    #colors = np.array([cmap(i) for i in np.linspace(first_rgb, last_rgb, ncolors).astype(int)])  # different depths of green
    colors = np.array([cmap(i) for i in [40, 110, 180]])  # different depths of green

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

    # Perform bootstrapping to estimate the the correlation coef and its uncertainty
    random.seed(2021)
    np.random.seed(2021)
    exp_h, exp_err = read_experimental_data()
    B25_avg, B25_std = sum_up_data(B25_sasa, B25_err)
    B26_avg, B26_std = sum_up_data(B26_sasa, B26_err)

    c1, _ = scipy.stats.kendalltau(B25_avg, exp_h)
    c2, _ = scipy.stats.kendalltau(B26_avg, exp_h)

    e1 = bootstrapping_sample(B25_sasa, exp_h, exp_err)
    e2 = bootstrapping_sample(B26_sasa, exp_h, exp_err)

    # Figure 1
    plt.figure(figsize=(6, 8))
    plt.subplot(2, 1, 1)
    for i in range(len(B25_avg)):
        if i == 0:
            plt.errorbar(B25_avg[i], exp_h[i], xerr=B25_std[i], yerr=exp_err[i], fmt="o", color='red', capsize=2, label='WT')
            plt.text(B25_avg[i] + 0.002, exp_h[i] + 0.1, sys[i], weight='bold')
        elif i == 1:  # just for the legend
            plt.errorbar(B25_avg[i], exp_h[i], xerr=B25_std[i], yerr=exp_err[i], fmt="o", color='blue', capsize=2, label='GF')
            plt.text(B25_avg[i] + 0.002, exp_h[i] + 0.1, sys[i])
        else:
            plt.errorbar(B25_avg[i], exp_h[i], xerr=B25_std[i], yerr=exp_err[i], fmt="o", color='blue', capsize=2)
            plt.text(B25_avg[i] + 0.002, exp_h[i] + 0.1, sys[i])
    
    plt.fill_between([0.02, 0.15], 2, 8, color=colors[0], alpha=0.3, zorder=0)
    plt.fill_between([0.02, 0.15], 8, 12.5, color=colors[1], alpha=0.3, zorder=0)
    plt.fill_between([0.02, 0.15], 12.5, 26, color=colors[2], alpha=0.3, zorder=0)
    plt.text(0.023, 12.8, '(Longer half-life than WT)')
    plt.text(0.023, 8.3, '(Comparable half-life as WT)')
    plt.text(0.023, 2.8, '(Shorter half-life than WT)')
    plt.text(0.76, 0.93, r'($\tau=$' + f'{c1:.3f} $\pm$ {e1:.3f})', transform=plt.gca().transAxes)
    
    plt.xlabel('SASA of B25-B26 scissile bond (nm$^2$)', size=12)
    plt.ylabel(r'$\alpha$-chymotrypsin half-life (min)', size=12)
    plt.xlim([0.02, 0.15])
    plt.ylim([2, 26])
    plt.grid()
    plt.legend(bbox_to_anchor=(0.98, 0.92))

    plt.text(-0.05, 1.20, 'A', transform=plt.gca().transAxes,
      fontsize=28, fontweight='bold', va='top', ha='right')

    # Figure 2
    plt.subplot(2, 1, 2)
    for i in range(len(B26_avg)):
        if i == 0:
            plt.errorbar(B26_avg[i], exp_h[i], xerr=B26_std[i], yerr=exp_err[i], fmt="o", color='red', capsize=2, label='WT')
            plt.text(B26_avg[i] + 0.002, exp_h[i] + 0.1, sys[i], weight='bold')
        elif i == 1:
            plt.errorbar(B26_avg[i], exp_h[i], xerr=B26_std[i], yerr=exp_err[i], fmt="o", color='blue', capsize=2, label='GF')
            plt.text(B26_avg[i] + 0.002, exp_h[i] + 0.1, sys[i])
        else:
            plt.errorbar(B26_avg[i], exp_h[i], xerr=B26_std[i], yerr=exp_err[i], fmt="o", color='blue', capsize=2)
            plt.text(B26_avg[i] + 0.002, exp_h[i] + 0.1, sys[i])
    
    plt.fill_between([0, 0.155], 2, 8, color=colors[0], alpha=0.3, zorder=0)
    plt.fill_between([0, 0.155], 8, 12.5, color=colors[1], alpha=0.3, zorder=0)
    plt.fill_between([0, 0.155], 12.5, 26, color=colors[2], alpha=0.3, zorder=0)
    plt.text(0.003, 12.8, '(Longer half-life than WT)')
    plt.text(0.003, 8.3, '(Comparable half-life as WT)')
    plt.text(0.003, 2.8, '(Shorter half-life than WT)')
    plt.text(0.76, 0.93, r'($\tau=$' + f'{c2:.3f} $\pm$ {e2:.3f})', transform=plt.gca().transAxes)
    
    plt.xlabel('SASA of B26-B27 scissile bond (nm$^2$)', size=12)
    plt.ylabel(r'$\alpha$-chymotrypsin half-life (min)', size=12)
    plt.xlim([0, 0.155])
    plt.ylim([2, 26])
    plt.grid()
    plt.legend(bbox_to_anchor=(0.98, 0.92))

    plt.tight_layout()
    #plt.suptitle(r'Correlation b/w $\alpha$-chymotrypsin half-life and scissile bond SASA', weight='bold', size=14)
    plt.savefig('sasa_PB_correlation.png', dpi=600)

    t2 = time.time()
    print(f"Elapsed time: {t2 - t1} seconds.")
    
