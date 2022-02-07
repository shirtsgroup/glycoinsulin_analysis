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

def locate_transition(xpm_data, r=500, N=10):
    """
    This function locates the rough time frame where the biggest transition occurs.
    """
    rmsd_1d = xpm_data[0]   # The rmsd value of the initial structure w.r.t all others
    diff = []
    for i in range(r, len(rmsd_1d) - r):
        avg_rmsd = [np.mean(rmsd_1d[i - r : i]), np.mean(rmsd_1d[i: i + r])]
        diff.append(max(avg_rmsd) / min(avg_rmsd)) 

    # Below we find the time frames corresopnding to the largest differences
    diff_idx = list(np.argsort(diff))  # idx of increasing values
    diff_idx.reverse()   # idx of decreasing values
    idx = np.array(diff_idx[:N])   # The indices of the N largest numbers

    transition = idx * (2000 / (len(xpm_data) - 1))  # ns

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
