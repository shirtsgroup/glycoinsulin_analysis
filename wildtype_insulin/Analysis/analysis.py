import os 

if __name__ == '__main__':
    # 1. RMSD calculaiton (select "Backbone" for least squares fit and RMSD calculation)
    os.system('echo 4 4 | gmx rms -s ../MD/4eyd_md.tpr -f ../MD/traj_comp.xtc -o 4eyd_rmsd_100ns.xvg')
    os.system('plot_2d -f 4eyd_rmsd_100ns.xvg -x "Time (ns)" -y "RMSD (nm)" -t "RMSD as a function of time" -n "4eyd_rmsd_100ns"')

    # 2. RMSF calculation (Group selected: "Backbone")
    os.system("echo 4 | gmx rmsf -f ../MD/traj_comp.xtc -s ../MD/4eyd_md.tpr -o 4eyd_rmsf_100ns.xvg -res")
    os.system('plot_2d -f 4eyd_rmsf_100ns.xvg -x "Residue" -y "RMSF (nm)" -t "RMSF as a function of time" -n "4eyd_rmsf_100ns"')

    # 3. Radus of gyration (Group selected: "Backbone")
    os.system('echo 4 | gmx gyrate -f ../MD/traj_comp.xtc -s ../MD/4eyd_md.tpr -o 4eyd_rGyr_100ns.xvg')
    os.system('plot_2d -f 4eyd_rGyr_100ns.xvg -x "Time (ns)" -y "Radius of gyration (nm)" -t "Radius or gyration as a function of time" -n "4eyd_rGyr_100ns"')

    # 4. SASA (Selected group for "surface": "Protein")
    os.system('echo 1 | gmx sasa -f ../MD/traj_comp.xtc -s ../MD/4eyd_md.tpr -o 4eyd_sasa_100ns.xvg')
    os.system('plot_2d -f 4eyd_sasa_100ns.xvg -x "Time (ns)" -y "Solvent accessible surface area ($ nm^{2} $)" -t "SASA as a function of time" -n "4eyd_sasa_100ns"')

    # 5. End-to-end distance (Selected group of polymer mainchain atoms: "Backbone")
    os.system('gmx make_ndx -f ../MD/4eyd_md.gro -o 4eyd.ndx')
    os.system('echo 19 | gmx polystat -s ../MD/4eyd_md.tpr -f ../MD/traj_comp.xtc -o 4eyd_eted_chainA_100ns.xvg -n 4eyd.ndx')
    os.system('echo 20 | gmx polystat -s ../MD/4eyd_md.tpr -f ../MD/traj_comp.xtc -o 4eyd_eted_chainB_100ns.xvg -n 4eyd.ndx')
    os.system('plot_2d -f 4eyd_eted_chainA_100ns.xvg -x "Time (ns)" -y "End-to-end distance (nm)" -t "End-to-end distance as a function of time (chain A)" -n "4eyd_eted_chainA_100ns"')
    os.system('plot_2d -f 4eyd_eted_chainA_100ns.xvg -x "Time (ns)" -y "End-to-end distance (nm)" -t "End-to-end distance as a function of time (chain B)" -n "4eyd_eted_chainB_100ns"')

    # 6. Number of hydrogen bonds (Two groups for analysis: "Protein" and "Protein")
    os.system('echo 1 1 | gmx hbond -f ../MD/traj_comp.xtc -s ../MD/4eyd_md.tpr -num 4eyd_hbond_100ns.xvg')
    os.system('plot_2d -f 4eyd_hbond_100ns.xvg -x "Time (ns)" -y "Numer of intramolecular hydrogon bonds" -t "Number of H.B. as a function of time" -n "4eyd_hbond_100ns"')

    # 7. DSSP (Selected group for analysis: "Protein")
    os.system('echo 1 | gmx do_dssp -f ../MD/traj_comp.xtc -s ../MD/4eyd_md.tpr -sc 4eyd_sscount_100ns.xvg -o 4eyd_dssp_100ns.xpm')
    os.system('')

    # 8. Animated trajectory 
    os.system("echo 1 0 | gmx trjconv -s ../MD/4eyd.tpr -f ../traj_comp.xtc -o 4eyd_nojump.xtc -center -pbc nojump")
    os.system("echo 1 0 | gmx trjconv -s ../MD/4eyd.tpr -f 4eyd_nojump.xtc -o 4eyd_center_dt250.xtc -center -pbc mol -ur compact -dt 250")


