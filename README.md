Glycoinsulin_project
=======================

## Description
This is a repository for maintaining the important files for running or analyzing MD simulations of different insulin wildtypes and their glycovariants. Specifically, we adopted five X-ray/NMR structures of the insulin wildtype, whose details for file preparation for MD simulations are summarized below. 

| PDB code      | 4EYD           | 4EY9           | 4EY1           | 3I3Z                   | 2MVC               |
|:-------------:|:--------------:|:--------------:|:--------------:|:----------------------:|:------------------:|
| PDB for H++   | 4eyd_clean.pdb | 4ey9_clean.pdb | 4ey1_clean.pdb | 3i3z_clean_mutated.pdb | 2mvc_clean_noH.pdb |
| Resolution    | 1.47 angstrom  | 1.47 angstrom  | 1.47 angstrom  | 1.60 angstrom          | 1.60 angstrom      |
| Method        | X-ray          | X-ray          | X-ray          | X-ray                  | NMR                |
| Preprocessing | Hpp            | Hpp            | Hpp            | PyMol + Hpp            | PyMol + Hpp        |
| pH value      | 8.0            | 8.0            | 7.9            | 6.9                    | 7.3                |
| Total charges | -1             | -1             | -1             | -1                     | -1                 |

Note that due to the huge size of the simulation files, here we only maintain the following files that are important for preparing/analyzing the MD simulations. We don't store the output files of the MD simulations. Specifically, in this repository, there are the following folders
- `wildtype_insulin`:
  - `Hpp_results/`: Output files of insulin wildtype obtained from the H++ server
- `glyco_insulin`:
  - `Glycam_outputs/`: Output files obtained from the GLYCAM glycoprotein builder.
- `useful_files`:
  - `acpype.py`: A Python code for converting coordinates and topology files between different the formats compatible with different simulation software. 
  -  `mdp_files`: mdp files for neutralization, energy minimization, equilibration and MD simulation. 
  - `prep_dir.sh`: A simple bash script for preparing the folders given the mdp files. 
  - `prep_md.sh`: A simple bash script for automating all the steps prior to the MD simulation. 
  - `gmx_analysis.sh`: A simple bash script calling GROMACS commands to perform the most commonly used data analysis.
- `analysis_codes`: 
  - Assessment of proteolytic stability
    - `SASA_analysis.py`: A Python code for calculating the solvent-accessible surface area (SASA) of the cleavage site for each glycoform. The results include histograms, time evolutions and sitributions of the SASA value of each glycoform. Note that the peptide bond plane (CONH) plane was considered when calculating the SASA, which involve the atoms of the cleavage site itself and its next residue. 
    - `Baker_Hubbard_hbonds.py`: A Python code for identifying the glycan-involved hydrogen bonds using the Baker-Hubbard criteria. The existence percentage of each hydrogen bond is calculated. 
  - Assessment of dimerization propensity
  - Assessment of biological activity


For more details about the project overview, protocols for file preparations, simulations and data analysis, please refer to our [research note hosted on HackMD](https://hackmd.io/@WeiTseHsu/glycoengineering).




