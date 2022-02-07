import numpy as np 
import gromacs as gmx 

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
    data = gmx.fileformats.xpm.XPM('4ey9_pairwise_rmsd.xpm')
