#!/bin/sh
set -e   # exit upon error 

echo The script launches a bunch of GROMACS command the perform the most commonly used data analysis for an MD simulation.

read -p "Are you working with in insulin wildypte (1) or an insulin glycoform (2)?" type
if [ ${type} == "1" ]
then
    read -p "Plase input the system name: " sys
elif [ ${type} == "2" ]
then 
    read -p "Please input the number of the ACS glycoform: " n
    sys=glycoform_${n}_ACS
fi

gmx_ver="gmx_mpi"
len="2000"

# 1. RMSD
echo 4 4 | mpirun -np 1 ${gmx_ver} rms -s ../MD/${sys}_md.tpr -f ../MD/traj_comp.xtc -o ${sys}_rmsd_${len}ns_dt500.xvg -dt 500 

# 2. RMSF
mpirun -np 1 ${gmx_ver} rmsf -f ../MD/traj_comp.xtc -s ../MD/${sys}_md.tpr -o ${sys}_rmsf_${len}ns_dt500.xvg -res -dt 500

# 3. Radius of gyration
mpirun -np 1 ${gmx_ver} gyrate -f ../MD/traj_comp.xtc -s ../MD/${sys}_md.tpr -o ${sys}_rGyr_${len}ns_dt500.xvg -dt 500

# 4. Animated trajectory
${gmx_ver} trjconv -s ../MD/${sys}_md.tpr -f ../MD/traj_comp.xtc -o ${sys}_${len}ns_nojump.xtc -center -pbc nojump
${gmx_ver} trjconv -s ../MD/${sys}_md.tpr -f ${sys}_${len}_nojump.xtc -o ${sys}_${len}ns_center_dt500.xtc -center -pbc mol -ur compact -dt 500

# Some additional analysis
# 5. End-to-end distance of chain A
#${gmx_ver} polystat -s ../MD/${sys}_md.tpr -f ../MD/traj_comp.xtc -o ${sys}_eted_chainA_${len}ns.xvg -n ${sys}.ndx

# 6. End-to-end distance of chain B
#${gmx_ver} polystat -s ../MD/${sys}_md.tpr -f ../MD/traj_comp.xtc -o ${sys}_eted_chainB_${len}ns.xvg -n ${sys}.ndx

# 7. Number of hydrogen bonds
#${gmx_ver} hbond -f ../MD/traj_comp.xtc -s ../MD/${sys}_md.tpr -num ${sys}_hbond_${len}ns.xvg

# 7. DSSP
#${gmx_ver} do_dssp -f ../MD/traj_comp.xtc -s ../MD/${sys}_md.tpr -sc ${sys}_sscount_${len}ns.xvg -o ${sys}_dssp_${len}ns.xpm

