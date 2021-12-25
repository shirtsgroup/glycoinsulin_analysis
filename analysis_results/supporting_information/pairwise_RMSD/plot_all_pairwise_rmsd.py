import time
import gromacs as gmx 
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc 
import matplotlib.gridspec as gridspec

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
