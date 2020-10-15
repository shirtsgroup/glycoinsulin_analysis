#!/bin/sh
#SBATCH --job-name insulin
#SBATCH -p RM
#SBATCH -N 4
#SBATCH -t 24:00:00
#SBATCH --ntasks-per-node=1

source /home/wehs7661/pkgs/gromacs/2020.2/bin/GMXRC
source /home/wehs7661/src/plumed2/sourceme.sh
module load gcc/10.1.0
module load mpi/intel_mpi

mpirun -np 16 gmx_mpi mdrun -s 4eyd_md.tpr -o 4eyd_md.trr -c 4eyd_md.gro -g md.log -e md.edr
