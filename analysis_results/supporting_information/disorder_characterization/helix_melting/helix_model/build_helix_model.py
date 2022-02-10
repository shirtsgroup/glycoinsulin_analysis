import os
import numpy as np
from cg_openmm.utilities.helix_optimize_nonbonded import *

# Optimize LJ sigma parameters of backbone and sidechain beads, for a fixed
# helical radius and pitch, and fixed bond lengths.

# Epsilon parameters are assumed equal for backbone and sidechain.

# References of the parameters:
# Ref 1: Characterization of Biomolecular Helices and Their Complementarity Using Geometric Analysis
# Ref 2: Minimialist models for proteins: A comparative analysis

# Helical parameters:
radius = 2.3 * unit.angstrom  # From Ref 1
pitch = 5.5 * unit.angstrom   # From Ref 1

# Bond lengths:
bond_dist_bb = 3.8 * unit.angstrom   # From Ref 2
bond_dist_sc = 3.8 * unit.angstrom   # This does not really matter in our case

# Number of backbone particles:
n_particle_bb = 9   # number of residues we want in the PDB file

opt_solution, geometry = optimize_helix_LJ_parameters(
    radius, pitch, n_particle_bb,
    bond_dist_bb=bond_dist_bb, bond_dist_sc=bond_dist_sc,
    pdbfile='helix_opt_LJ_params_constrained.pdb',
    plotfile='helix_opt_LJ_params_constrained.pdf',
    DE_popsize=50)

print(opt_solution)
print(geometry)
