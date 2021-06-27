#!/bin/sh
set -e   # exit upon error
echo This script automates all the steps before launch an MD simulation.

# Step 1: Take in the system name 
read -p "Are you working with in an insulin wildypte (1), an insulin glycoform in the 2018 ACS paper (2) or other (3)? " type
if [ ${type} == "1" ]
then
    read -p "Plase input the system name: " sys
    mkdir Box EM Equil MD Sol_ions Topology Analysis
    cd Topology
    cp ../Hpp_results/${sys}.top .
    cp ../Hpp_results/${sys}.gro .
elif [ ${type} == "2" ]
then 
    read -p "Please input the serial number of the ACS glycoform: " n
    sys=glycoform_${n}_ACS
    mkdir Box EM Equil MD Sol_ions Topology Analysis
    cd Topology
    cp ../Glycam_outputs/structure_GMX.gro ${sys}.gro
    cp ../Glycam_outputs/structure_GMX.top ${sys}.top
elif [ ${type} == "3" ]
then 
    read -p "Plase input the system name: " sys
    mkdir Box EM Equil MD Sol_ions Topology Analysis
    cd Topology
    cp ../Glycam_outputs/structure_GMX.gro ${sys}.gro
    cp ../Glycam_outputs/structure_GMX.top ${sys}.top
fi

# Step 2: Modify the topology and distribute the mdp files
d1=" OW       OW          16.00    0.0000    A     3.15061e-01  6.36386e-01
 HW       HW          1.008    0.0000    A     0.00000e+00  0.00000e+00
 Cl       Cl          35.45    0.00000   A    4.40104e-01    4.18400e-01
 Na       Na          22.99    0.00000   A    3.32840e-01    1.15897e-02
 IB       IB          0.00000  0.00000   A    8.90899e-01    4.18400e-01
 C0       C0          0.00000  0.00000   A    3.05240e-01    1.92376e+00
 MG       MG          0.00000  0.00000   A    1.41225e-01    3.74342e+00
 K        K           0.00000  0.00000   A    4.73602e-01    1.37235e-03
 Rb       Rb          0.00000  0.00000   A    5.26699e-01    7.11280e-04
 Cs       Cs          0.00000  0.00000   A    6.04920e-01    3.37230e-04
 Li       Li          0.00000  0.00000   A    2.02590e-01    7.65672e-02
 Zn       Zn          0.00000  0.00000   A    1.95998e-01    5.23000e-02"

d2="; Include water topology
#include \"amber99sb-ildn.ff/tip3p.itp\"

; Include topology for ions
#include \"amber99sb-ildn.ff/ions.itp\"
"

preprocessed_d1=$(printf '%s\n' "$d1" |
  sed 's/\\/&&/g;s/^[[:blank:]]/\\&/;s/$/\\/')

preprocessed_d2=$(printf '%s\n' "$d2" |
  sed 's/\\/&&/g;s/^[[:blank:]]/\\&/;s/$/\\/')

sed "/;name   bond_type     mass     charge   ptype   sigma         epsilon       Amb/a ${preprocessed_d1%?}" ${sys}.top -i.bkp1
sed "/\[ system \]/i ${preprocessed_d2%?}" ${sys}.top -i.bkp2

cd ../ && mv gmx_analysis.sh Analysis
cd Sol_ions && cp ../mdp_files/ions.mdp .
cd ../EM && cp ../mdp_files/em.mdp .
cd ../Equil && mkdir NVT NPT 
cd NVT && cp ../../mdp_files/nvt_equil.mdp .
cd ../NPT && cp ../../mdp_files/npt_equil.mdp .
cd ../../MD && cp ../mdp_files/md.mdp .


# Step 3: Run all the steps before launching the MD simulation
cd ../Box && mpirun -np 1 gmx_mpi editconf -f ../Topology/${sys}.gro -o ${sys}_box.gro -bt dodecahedron -d 1.0 -c

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

echo Given the salinity as 0.15 M, below is the number of ions required given the box volume:
table="+-----------+---------------------+
| # of ions |      Box volume     |
+-----------+---------------------+
|     10    |  99.6323 - 110.7025 |
|     11    | 110.7026 - 121.7728 |
|     12    | 121.7729 - 132.8431 |
|     13    | 132.8432 - 143.9133 |
|     14    | 143.9134 - 154.9836 |
+-----------+---------------------+"
echo "${table}"
qtot=$(grep "qtot" ../Topology/${sys}.top | tail -1 | awk '{print $11}')
echo "According to the topology file, the number of total charges is ${qtot}"
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

read -p "The job name is the same as the system name (${sys}). Do you wnat to change it? If so, enter a new name, or, enter 'no': " change
if [ ${change} == "no" ]
then
    jobname=${sys}
else
    jobname=${change}
fi

echo Writing job submission script ...
touch run.sh
d="#!/bin/sh
#SBATCH --job-name ${jobname}
#SBATCH -p RM
#SBATCH -N 1
#SBATCH -t 48:00:00
#SBATCH --ntasks-per-node=128

source /jet/home/wehs7661/pkgs/gromacs/2020.4/bin/GMXRC
source /jet/home/wehs7661/src/plumed2/sourceme.sh
module load gcc/10.2.0
module load openmpi/3.1.6-gcc10.2.0

mpirun -np 128 gmx_mpi mdrun -s ${sys}_md.tpr -x ${sys}_md.xtc -c ${sys}_md.gro -g md.log -e md.edr -cpi state.cpt -ntomp 1 -npme 32"
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
