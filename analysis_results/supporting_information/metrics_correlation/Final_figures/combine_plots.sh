#combine_plots -f ../figures/sasa_sc_B25__sasa_res_B24.png ../figures/sasa_sc_B25__sasa_res_B25.png ../figures/sasa_sc_B25__sasa_sc_B26.png ../figures/sasa_sc_B26__sasa_res_B24.png ../figures/sasa_sc_B26__sasa_res_B25.png ../figures/sasa_res_B24__sasa_res_B25.png -n Fig_SASA_metrics_correlation -s 16 8

#combine_plots -f ../figures/beta_B22__beta_B23.png ../figures/beta_B22__beta_B24.png ../figures/beta_B22__beta_B25.png ../figures/beta_B23__beta_B24.png ../figures/beta_B23__beta_B25.png ../figures/beta_B24__beta_B25.png -n Fig_beta_metrics_correlation -s 16 8

combine_plots -f ../figures/mixed_metrics/* -n Fig_SASA_beta_metrics_correlation -s 20 16
