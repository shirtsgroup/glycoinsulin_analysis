#!/bin/sh
set -e   # exit upon error
echo This script create the folders and parameter files for later GROMACS commands required for all the steps before running the MD simulation.

read -p "Are you working with in insulin wildypte (1) or an insulin glycoform (2)?" type
if [ ${type} == "1" ]
then
    read -p "Plase input the system name: " sys
    mkdir Box EM Equil MD Sol_ions Topology Analysis
    cd Topology
    cp ../Hpp_results/${sys}.top .
    cp ../Hpp_results/${sys}.gro .
elif [ ${type} == "2" ]
then 
    read -p "Please input the number of the ACS glycoform: " n
    sys=glycoform_${n}_ACS
    mkdir Box EM Equil MD Sol_ions Topology Analysis
    cd Topology
    cp ../Glycam_outputs/structure_GMX.gro ${sys}.gro
    cp ../Glycam_outputs/structure_GMX.top ${sys}.top
fi

mv gmx_analysis.sh Analysis
cd ../Sol_ions && cp ../mdp_files/ions.mdp .
cd ../EM && cp ../mdp_files/em.mdp .
cd ../Equil && mkdir NVT NPT 
cd NVT && cp ../../mdp_files/nvt_equil.mdp .
cd ../NPT && cp ../../mdp_files/npt_equil.mdp .
cd ../../MD && cp ../mdp_files/md.mdp .


