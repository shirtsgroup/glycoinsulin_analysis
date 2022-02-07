import numpy as np 
import gromacs as gmx 

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
    data = gmx.fileformats.xpm.XPM('3i3z_pairwise_rmsd.xpm')
    rmsd = data.array.astype('float64')
    t = locate_transition(rmsd)
    print(t)

