# A quick guide for generating the figures in the paper
## 1. Figures in the main text
- Figure 1 to Figure 3 were made by PyMol and Illustrator. 
- Figure 4 (Correlation plot between SASA of scissile bonds/P1 sites and chymotrypsin half-life)
  - Figure 4A
    - Working directory: `Glycoinsulin_project/analysis_results/proteolytic_stability/Metric_1`
    - Command: `python plot_PB_SASA_correlation.py`
  - Figure 4B 
    - Working directory: `Glycoinsulin_project/analysis_results/proteolytic_stability/Metric_2`
    - Command: `python plot_res_SASA_correlation.py`
  - Figure 4
    - Working directory: `Glycoinsulin_project/analysis_results/proteolytic_stability/Final_figures`
    - Command: `combine_plots -f ../Metric_1/sasa_PB_correlation.png ../Metric_2/sasa_res_correlation.png -n Fig_SASA_correlation -s 15 10`
- Figure 5 (Correlation plot between beta-sheet propensity and chymotrypsin half-life)
  - Working directory: `Glycoinsulin_project/analysis_results/proteolytic_stability/Metric_3`
  - Command: `python plot_avg_beta_propensity_correlation.py`
- Figure 6 (Glycan-involved hydrogen bond)
  - Working directory: `Glycoinsulin_project/analysis_results/proteolytic_stability/Metric_4`
  - Command: `python hbond_plot.py`
- Figure 7 (Glycan-dimer occlusion)

## 2. Figures in the SI