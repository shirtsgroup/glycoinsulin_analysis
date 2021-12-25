File name: DSSP_averages_total_Fig7.py
	Description: This file has all the data necessary for plotting the DSSP data for the main text. Simply run the script and go. Do not alter the data within this file.
	This is the only file necessary for reproducing the graph in Figure 7 AND the graph used in Supplemental Figure S5.
	
File name: plotting_dsspBinomials.py
	Description: This file plots the individual model DSSP binomial data for the supplement. Two additional files are required, plus an additional "prefix" name.
	Usage: python3 plotting_dsspBinomials.py {model_totals.csv} {model_lagtimes.csv} {model_name}
	This is the script that makes individual graphs for Supplemental Figure S4. Each glycoform will need to be run individually to make all the graphs for the supplemental figure.
	
File name: suppFigS3_fractionDSSP_plotting.py
	Description: This file plots the fraction of total secondary structure calculated by DSSP for the supplemental figure S3. One additional file is required.
	Usage: python3 suppFigS3_fractionDSSP_plotting.py {model_totals.csv} {model_name}
		e.g. python3 suppFigS3_fractionDSSP_plotting.py 4eyd_totals.csv 4EYD
	This is the script necessary to make the individual graphs for Supplemental Figure S3. Each glycoform will need to be run to make all the graphs for this supplemental figure.
	
File name: suppFigS3_WT_fractionDSSP_plotting.py
	Description: This file plots the wild type graph data for Supplemental Figure 3. It does not require any additional files to run (it has the data hard coded).
	Usage: python3 suppFigS3_WT_fractionDSSP_plotting.py
	