import os 
import numpy as np
import mdtraj as md 
from prettytable import PrettyTable

def logger(*args, **kwargs):
    print(*args, **kwargs)
    with open('characterize_melting_results.txt', "a") as f:
        print(file=f, *args, **kwargs)

if __name__ == "__main__":
    # Load in and slice the helix model and the trajectories
    traj_tot = md.load('all_traj_CA.xtc', top='all_traj_CA.pdb')
    ref = md.load('helix_model/helix_model.pdb')
    s = int(traj_tot.n_frames / 5)   # size of each trajectory

    # Get the 4 segments of alpha helix (A1-A6, A2-A7, A3-A8, A4-A9)
    trajs = [traj_tot.atom_slice(range(i, i + 6)) for i in range(3)]
    refs = [ref.atom_slice(range(i, i + 6)) for i in range(3)]
    
    # Caculate H_segment and the fraction of helix
    H_seg = 0
    for i in range(3):
        rmsd = md.rmsd(trajs[i], refs[i])
        H_seg += (1 - np.power(rmsd, 8)) / (1 - np.power(rmsd, 12))
    
    x = PrettyTable()
    x.field_names = ['WT model', 'AN-helix melting']
    sys = ['4EYD', '4EY9', '4EY1', '3I3Z', '2MVC']
    for i in range(5):
        H = H_seg[s * i : s * (i + 1)]
        fraction = sum(H < 2) / len(H) * 100
        results = [f'{sys[i]}', f'{fraction:.2f}%']
        x.add_row(results)

    f = sum(H_seg < 2) / len(H_seg) * 100
    r = ['Total', f'{f:.2f}%']
    x.add_row(r)

    logger(x)
