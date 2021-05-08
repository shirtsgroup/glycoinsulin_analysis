set -e 
sys=4eyd
gmx_ver="gmx_mpi"

# 1. RMSD
mpirun -np 1 ${gmx_ver} rms -s ../MD/${sys}_md.tpr -f ../MD/traj_comp.xtc -o ${sys}_rmsd_500ns_dt500.xvg -dt 500 -e 500000
mpirun -np 1 ${gmx_ver} rms -s ../MD/${sys}_md.tpr -f ../MD/traj_comp.xtc -o ${sys}_rmsd_2000ns_dt500.xvg -dt 500 

#ï2. RMSF
mpirun -np 1 ${gmx_ver} rmsf -f ../MD/traj_comp.xtc -s ../MD/${sys}_md.tpr -o ${sys}_rmsf_500ns_dt500.xvg -res -dt 500 -e 500000
mpirun -np 1 ${gmx_ver} rmsf -f ../MD/traj_comp.xtc -s ../MD/${sys}_md.tpr -o ${sys}_rmsf_2000ns_dt500.xvg -res -dt 500

# 3. Radius of gyration
mpirun -np 1 ${gmx_ver} gyrate -f ../MD/traj_comp.xtc -s ../MD/${sys}_md.tpr -o ${sys}_rGyr_500ns_dt500.xvg -dt 500 -e 500000
mpirun -np 1 ${gmx_ver} gyrate -f ../MD/traj_comp.xtc -s ../MD/${sys}_md.tpr -o ${sys}_rGyr_2000ns_dt500.xvg -dt 500

# 4. End-to-end distance of chain A
#${gmx_ver} polystat -s ../MD/${sys}_md.tpr -f ../MD/traj_comp.xtc -o ${sys}_eted_chainA_${len}.xvg -n ${sys}.ndx

# 5. End-to-end distance of chain B
#${gmx_ver} polystat -s ../MD/${sys}_md.tpr -f ../MD/traj_comp.xtc -o ${sys}_eted_chainB_${len}.xvg -n ${sys}.ndx

# 6. Number of hydrogen bonds
#${gmx_ver} hbond -f ../MD/traj_comp.xtc -s ../MD/${sys}_md.tpr -num ${sys}_hbond_${len}.xvg

# 7. DSSP
#${gmx_ver} do_dssp -f ../MD/traj_comp.xtc -s ../MD/${sys}_md.tpr -sc ${sys}_sscount_100ns.xvg -o ${sys}_dssp_${len}.xpm

# 8. Animated trajectory
#${gmx_ver} trjconv -s ../MD/${sys}_md.tpr -f ../MD/traj_comp.xtc -o ${sys}_${len}_nojump.xtc -center -pbc nojump
#${gmx_ver} trjconv -s ../MD/${sys}_md.tpr -f ${sys}_${len}_nojump.xtc -o ${sys}_${len}_center_dt250.xtc -center -pbc mol -ur compact -dt 250

