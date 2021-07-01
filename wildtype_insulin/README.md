# Wildtype insulin structures
In the 5 subfolders here, including `2MVC`, `3I3Z`, `4EY1`, `4EY9`, and `4EYD`, we maintain the input (`{PDB_ID}_clean*.pdb`) and output files of the H++ server (all stored in the folder `Hpp_results`) for preparing different structures of the insulin wildtype. Note that each of the wildtype structures were prepared at different pH values to have the same amount of total charges. Below we summarize some details of preparing these structures for the MD simulations. 

| PDB code      | 4EYD           | 4EY9           | 4EY1           | 3I3Z                   | 2MVC               |
|:-------------:|:--------------:|:--------------:|:--------------:|:----------------------:|:------------------:|
| PDB for H++   | 4eyd_clean.pdb | 4ey9_clean.pdb | 4ey1_clean.pdb | 3i3z_clean_mutated.pdb | 2mvc_clean_noH.pdb |
| Resolution    | 1.47 angstrom  | 1.47 angstrom  | 1.47 angstrom  | 1.60 angstrom          | 1.60 angstrom      |
| Method        | X-ray          | X-ray          | X-ray          | X-ray                  | NMR                |
| Preprocessing | Hpp            | Hpp            | Hpp            | PyMol + Hpp            | PyMol + Hpp        |
| pH value      | 8.0            | 8.0            | 7.9            | 6.9                    | 7.3                |
| Total charges | -1             | -1             | -1             | -1                     | -1                 |

