import time
import gromacs as gmx 
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc 
import matplotlib.gridspec as gridspec

def logger(*args, **kwargs):
    print(*args, **kwargs)
    with open('transition_frames.txt', "a") as f:
        print(file=f, *args, **kwargs)

def locate_transition(xpm_data, w):
    """
    This function locates the rough time frame where the biggest transition occurs.
    For more information about this method, check change_point_detection.ipynb
    """
    rmsd_1d = xpm_data[0]   # The rmsd value of the initial structure w.r.t all others
    std_ratio = []
    for i in range(w, len(rmsd_1d) - w):  # w: window size
        window_std = [np.std(rmsd_1d[i - w : i]), np.std(rmsd_1d[i : i + w])]
        std_ratio.append(max(window_std) / min(window_std))
    
    idx = std_ratio.index(max(std_ratio)) + w
    transition = idx * (2000 / (len(xpm_data) - 1))  # ns

    return transition, std_ratio

if __name__ == "__main__":
    t1 = time.time()

    rc('font', **{
       'family': 'sans-serif',
       'sans-serif': ['DejaVu Sans'],
       'size': 10
    })
    # Set the font used for MathJax - more on this later
    rc('mathtext', **{'default': 'regular'})
    plt.rc('font', family='Arial')

    w_size = [100, ]  # values from the preliminary tests in change_point_detection.ipynb
    sys = ['4eyd', '4ey9', '4ey1', '3i3z', '2mvc']
    labels = ['A', 'B', 'C', 'D', 'E']
    fig_name = 'all_pairwise_rmsd.png'

    fig = plt.figure(figsize=(12, 7), tight_layout=True)
    spec = gridspec.GridSpec(ncols=6, nrows=2, figure=fig)
    specs = [spec[0, 0:2], spec[0, 2:4], spec[0, 4:6], spec[1, 1:3], spec[1, 3:5]]
    for i in range(len(sys)):
        xpm = f'{sys[i]}_pairwise_rmsd.xpm'
        data = gmx.fileformats.xpm.XPM(xpm)
        rmsd = data.array.astype('float64')

        # Identify the time frame where the transition occurred
        w = 100
        t, std_ratio = locate_transition(rmsd, w=w_size[i])
        logger(f'{sys[i].upper()}: The largest transition between states occurred at {t} ns.')

        # plot RMSD w.r.t the initial structure
        plt.figure(figsize=(8, 6))
        plt.subplot(2, 1, 1)
        plt.plot(np.arange(len(rmsd)) * 0.25, rmsd[0])
        plt.xlabel('The time frame of the reference (ns)')
        plt.ylabel('RMSD w.r.t. the initial structure')
        plt.xlim([0, 2000])
        plt.grid()

        plt.subplot(2, 1, 2)
        plt.plot(np.arange(w, len(rmsd) - w) * 0.25, std_ratio)
        plt.xlabel('Time (ns)')
        plt.ylabel('s.t.d. ratio')
        plt.xlim([0, 2000])
        plt.grid()

        plt.tight_layout()
        plt.savefig(f'RMSD_initial_structure/{sys[i]}_rmsd_initial_structure.png', dpi=600)

        # plot pairwise RMSD
        ax = fig.add_subplot(specs[i])
        plt.imshow(rmsd, cmap='YlGnBu', interpolation="none", extent=[0, 2000, 2000, 0])
        plt.xlabel('Time (ns)')
        plt.ylabel('Time (ns)')
        plt.title(f'Pairwise RMSD of {sys[i].upper()} wildtype', weight='bold')
        plt.text(-0.28, 1.05, f'{labels[i]}', transform=plt.gca().transAxes, weight='bold', size=15)
        plt.colorbar(label='RMSD (nm)', fraction=0.05, pad=0.04)
        plt.gca().invert_yaxis()
    plt.subplots_adjust(wspace=0, hspace=0)
    #plt.tight_layout()

    plt.savefig(fig_name, dpi=300)

    t2 = time.time()
    print(f'Elapsed time: {t2 - t1} seconds.')
