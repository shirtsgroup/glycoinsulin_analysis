#!/bin/sh
set -e   # exit upon error
echo This script automates all the steps required before running a vanill MD simulation, including the construction of the box, solvation and neutralization of the system, energy minimization and equilibrations. This script also submit the job to run the MD simulation if needed. 

read -p "Are you working with in insulin wildypte (1) or an insulin glycoform (2)?" type
if [ ${type} == "1" ]
then
    read -p "Plase input the system name: " sys
elif [ ${type} == "2" ]
then 
    read -p "Please input the number of the ACS glycoform: " n
    sys=glycoform_${n}_ACS

cd Box && mpirun -np 1 gmx_mpi editconf -f ../Topology/${sys}.gro -o ${sys}_box.gro -bt dodecahedron -d 1.0 -c

cd ../Sol_ions && mpirun -np 1 gmx_mpi solvate -cp ../Box/${sys}_box.gro -p ../Topology/${sys}.top -o ${sys}_sol.gro -cs
mpirun -np 1 gmx_mpi grompp -f ions.mdp -c ${sys}_sol.gro -p ../Topology/${sys}.top -o ${sys}_ions.tpr -maxwarn 1

read -p "Do you want to stop and tweak the charges (yes/no): " stop
if [ ${stop} == "yes" ]
then
    exit
elif [ ${stop} == "no" ]
then 
    :
fi

echo Given the salinity as 0.15 M, as long as the box volume is between 110.7026 and 121.7728 nm^3, 11 ions from the buffer solution should be added.
read -p "Please input the number of sodium ions: " np
read -p "Please input the number of chloride ions: " nn
mpirun -np 1 gmx_mpi genion -s *tpr -o ${sys}_ions.gro -p ../Topology/${sys}.top -pname NA -nname CL -np ${np} -nn ${nn}
cd ../EM && mpirun -np 1 gmx_mpi grompp -f em.mdp -c ../Sol_ions/${sys}_ions.gro -p ../Topology/${sys}.top -o ${sys}_em.tpr -maxwarn 1
mpirun -np 64 gmx_mpi mdrun -s *tpr -o ${sys}_em.trr -c ${sys}_em.gro -g em.log -e em.edr -ntomp 1
cd ../Equil/NVT && mpirun -np 1 gmx_mpi grompp -f nvt_equil.mdp -c ../../EM/${sys}_em.gro -p ../../Topology/${sys}.top -o ${sys}_equil.tpr
mpirun -np 64 gmx_mpi mdrun -s *tpr -o ${sys}_equil.trr -c ${sys}_equil.gro -g equil.log -e equil.edr -ntomp 1
cd ../NPT/ && mpirun -np 1 gmx_mpi grompp -f npt_equil.mdp -c ../NVT/${sys}_equil.gro -t ../NVT/state.cpt -p ../../Topology/${sys}.top -o ${sys}_equil.tpr
mpirun -np 64 gmx_mpi mdrun -s *tpr -o ${sys}_equil.trr -c ${sys}_equil.gro -g equil.log -e equil.edr -ntomp 1
cd ../../MD/ && mpirun -np 1 gmx_mpi grompp -f md.mdp -c ../Equil/NPT/${sys}_equil.gro -p ../Topology/${sys}.top -o ${sys}_md.tpr -t ../Equil/NPT/state.cpt

echo Writing job submission script and submitting the job ...
touch run.sh
d="#!/bin/sh
#SBATCH --job-name ${sys}
#SBATCH -p RM
#SBATCH -N 1
#SBATCH -t 48:00:00
#SBATCH --ntasks-per-node=128

source /jet/home/wehs7661/pkgs/gromacs/2020.4/bin/GMXRC
source /jet/home/wehs7661/src/plumed2/sourceme.sh
module load gcc/10.2.0
module load openmpi/3.1.6-gcc10.2.0

mpirun -np 128 gmx_mpi mdrun -s ${sys}_md.tpr -o ${sys}_md.trr -c ${sys}_md.gro -g md.log -e md.edr -ntomp 1 -npme 32"
echo "${d}" >> run.sh

read -p "Do you want to submit the job (yes/no): " submit
if [ ${submit} == "yes" ]
then
    sbatch run.sh
    echo Job submitted!
elif [ ${submit} == "no" ]
then 
    exit
fi


