## Description
This folder stores data analysis scripts and results.
- `archved_scripts`: A folder containing deprecated scripts for data analysis.
- `dimerization_propensity`: A folder containing scripts for assessing the dimerization propensity of glyco-insulin and the corresponding results. 
  - `Metric_1`: Glycan-dimer occlusion analysis 
  - `Metric_1_deprecated`: DSSP analysis. Deprecated and not presented in the final work.
- `proteolytic_stability`:
  - `Exp_data`: Experimental data of alpha-chymotrypsin half-life of each glyco-variant.
  - `Final_figures`: Figures reported in the main text.
  - `Metric_1`: Calculation of scissile bond SASA analysis
  - `Metric_2`: Calculation of P1-site SASA
  - `Metric_3`: Calculation of beta-sheet propensity of the P1-P3 region
  - `Metric_4`: Glycan-involved hydrogen bond analysis 
- `raw_data_generation`: Python scripts for generating GROMACS `xvg` files from MD trajectories. 
- `supporting_information`: A folder containing scripts for performing additional analysis and the corresponding results, which has the following subfolders:
  - `disorder_characterization`: Characterization of three common disorder elements, including AN-helix melting, detachment of B-chain N-temrinus and detachment of B-chain C-terminus
  - `metric_correlation`: Correlation between 8 measures in Metric 1 to Metric 3 for the proteolytic stability
  - `pairwise RMSD`: Calculation of pairwise RMSD to demonstrate the sampling of the configurational space
  - `SASA_distributions`: KDE distribution of Metric 1 and Metric 2 of the proteolytic stability
  - `transition analysis`: Repeated data analysis for the states before and after transition in the pairwise RMSD values.
