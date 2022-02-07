import time
import argparse
import gromacs as gmx 
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib import rc 

def initialize():
    parser = argparse.ArgumentParser(
        description='This code reads in the xpm file containing pairwise rmsd data and visualize it.')
    parser.add_argument('-x',
                        '--xpm',
                        help='The filename of the xpm file.')
    parser.add_argument('-t',
                        '--time',
                        help='The time frame to use as the reference to output the 1D RMSD time series.')
    args_parse = parser.parse_args()

    return args_parse

def locate_transition(xpm_data, r=10):
    """
    This function locates the rough time frame where the biggest transition occurs.
    """
    rmsd_diff = []
    idx = range(r, len(xpm_data) - r)
    for i in idx:
        # In the case that range=10, for time frame i, we average the RMSD of the corresponding structure w.r.t. the 
        # structure at i-1, i-2, ..., i-10 (rsmd_1) and at i + 1, i + 2, ..., i + 10 (rmsd_2)
        # The time frame having the biggest difference between rmsd_1 and rmsd_2 is what we look for. 
        rmsd_1 = np.mean(xpm_data[i - r : i - 1, i])
        rmsd_2 = np.mean(xpm_data[i + 1 : i + r, i])
        rmsd_diff.append(np.abs(rmsd_1 - rmsd_2))

    transition_idx = idx[rmsd_diff.index(np.max(rmsd_diff))]
    print(rmsd_diff[transition_idx])
    transition = transition_idx * (2000 / (len(xpm_data) - 1))   # ns
    
    return transition


if __name__ == "__main__":
    t1 = time.time()
    args = initialize()

    rc('font', **{
       'family': 'sans-serif',
       'sans-serif': ['DejaVu Sans'],
       'size': 10
    })
    # Set the font used for MathJax - more on this later
    rc('mathtext', **{'default': 'regular'})
    plt.rc('font', family='Arial')

    sys = args.xpm.split('_')[0].upper()
    fig_name = args.xpm.split('.xpm')[0] + '.png'

    data = gmx.fileformats.xpm.XPM(args.xpm)
    rmsd = data.array.astype('float64')
    plt.figure()
    #sns.heatmap(rmsd, cmap="YlGnBu", cbar_kws={'label': 'RMSD (nm)'})
    plt.imshow(rmsd, cmap='YlGnBu', interpolation="none", extent=[0, 2000, 2000, 0])
    plt.xlabel('Time (ns)')
    plt.ylabel('Time (ns)')
    plt.title(f'Pairwise RMSD of {sys} wildtype', weight='bold')
    plt.colorbar(label='RMSD (nm)')
    plt.gca().invert_yaxis()
    plt.tight_layout(rect=[-0.2, 0, 1, 1])
    plt.savefig(fig_name, dpi=600)

    t2 = time.time()
    print(f'Elapsed time: {t2 - t1} seconds.')
