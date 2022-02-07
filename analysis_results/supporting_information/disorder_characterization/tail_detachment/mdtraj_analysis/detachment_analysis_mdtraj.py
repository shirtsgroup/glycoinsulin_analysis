import os
import time
import datetime
import numpy as np
import mdtraj as md 
from prettytable import PrettyTable

def logger(*args, **kwargs):
    print(*args, **kwargs)
    with open('detachment_analysis_results.txt', "a") as f:
        print(file=f, *args, **kwargs)

def format_time(t):
    hh_mm_ss = str(datetime.timedelta(seconds=t)).split(':')
    hh, mm, ss = float(hh_mm_ss[0]), float(hh_mm_ss[1]), float(hh_mm_ss[2])
    if hh == 0:
        if mm == 0:
            t_str = f'{ss} second(s)'
        else:
            t_str = f'{mm:.0f} minute(s) {ss:.0f} second(s)'
    else:
        t_str = f'{hh:.0f} hour(s) {mm:.0f} minute(s) {ss:.0f} second(s)'

    return t_str

def get_disorder_fraction(theta_1, theta_2, theta_3):
    # This funciton returns the percentage of each disorder element
    n_b1b7_flip = sum(theta_2 > 10) / len(theta_2) * 100
    n_b1b7_detach = sum(theta_1[theta_2 < 10] > 85) / len(theta_1) * 100
    n_b20b30_detach = sum(theta_3 < 0) / len(theta_3) * 100

    return n_b1b7_flip, n_b1b7_detach, n_b20b30_detach

if __name__ == "__main__":
    t1 = time.time()
    # Working directory: /ocean/projects/cts160011p/wehs7661/Glycoinsulin_project/wildtype_insulin

    # Step 1: Load in the trajectories
    sys = ['4eyd', '4ey9', '4ey1', '3i3z', '2mvc']
    if os.path.isfile(f'all_traj_CA.xtc') is True:
        traj_tot = md.load('all_traj_CA.xtc', top='all_traj_CA.pdb')
    else:
        paths = []
        prjt_dir = '/ocean/projects/cts160011p/wehs7661/Glycoinsulin_project/wildtype_insulin/'
        paths.append(prjt_dir + '4EYD/pH_8.0/MD/')
        paths.append(prjt_dir + '4EY9/pH_8.0/MD/')
        paths.append(prjt_dir + '4EY1/pH_7.9/MD/')
        paths.append(prjt_dir + '3I3Z/pH_6.9/MD/')
        paths.append(prjt_dir + '2MVC/pH_7.3/MD/')

        trajs = []
        for i in range(len(paths)):
            logger(f'Loading in to trajectory of the wild-type model {sys[i].upper()} ...')
            traj = md.load(paths[i] + f'{sys[i]}_md.xtc', top = paths[i] + f'{sys[i]}_md.pdb', stride=100)  # 10001 frames
            CA_idx = traj.top.select('name CA')
            traj = traj.atom_slice(CA_idx)
            trajs.append(traj) 
        traj_tot = md.join(trajs)
        
        traj.save('all_traj_CA.xtc')
        pdb = md.load(paths[0] + f'{sys[0]}_md.pdb')
        pdb = pdb.atom_slice(CA_idx)
        pdb.save('all_traj_CA.pdb')

    # Step 2: Calculate relevant angles
    res = [[23, 29, 40], [23, 14, 38, 35], [12, 18, 39, 45]]  # for B1-B7 detachment, B1-B7 flipping and B20-B30 detachment
    
    # 2-1: angle b/w CA atoms of B3, B9, and B20, for B1-B7 detachment
    angle = md.compute_angles(traj_tot, np.array([res[0]])) * 180 / np.pi

    # 2-2: dihedrals for B1-B7 flipping and B20-B30 detachment
    dihedral = md.compute_dihedrals(traj_tot, np.array(res[-2:])) * 180 / np.pi

    # Step 3: Calculate the fraction of each disorder element
    n_frames = traj_tot.n_frames
    theta_1 = angle.T.reshape(n_frames,)     # B1-B7 alignment: theta_1 > 85 and theta_2 < 10
    theta_2 = dihedral.T[0]               # B1-B7 flipping: theta_2 > 10
    theta_3 = dihedral.T[1]               # B20-B30 alignment: theta_3 < 0

    s = int(n_frames / len(sys))  # size of each trajectory
    x = PrettyTable()
    x.field_names = ['WT model', 'B1-B7 flipping', 'B1-B7 detachment', 'B20-B30 detachment']
    for i in range(5):
        fractions = get_disorder_fraction(theta_1[s * i : s * (i + 1)], theta_2[s * i : s * (i + 1)], theta_3[s * i : s * (i + 1)])
        results = [f'{sys[i]}'.upper()]
        results.extend([f'{fractions[i]:.2f}%' for i in range(len(fractions))])
        x.add_row(results)
    
    f = get_disorder_fraction(theta_1, theta_2, theta_3)
    r = ['Total']
    r.extend(f'{f[i]:.2f}%' for i in range(len(f)))
    x.add_row(r)

    logger(x)

    t2 = time.time()
    logger(f'Elapsed time: {format_time(t2 - t1)} minutes')