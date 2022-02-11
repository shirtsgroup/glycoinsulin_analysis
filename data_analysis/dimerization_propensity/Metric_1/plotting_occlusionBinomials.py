#!/usr/bin/env python3

import matplotlib.pyplot as plt 
from matplotlib import rc
import numpy as np
import argparse
import math
import pandas as pd

def wilsonScore(p, n, z):
    score_plus = (p + (z**2/(2*n)) + z*math.sqrt((p*(1-p)/n) + (z**2/(4*n**2)))) / (1 + (z**2/n))
    score_minus = (p + (z**2/(2*n)) - z*math.sqrt((p*(1-p)/n) + (z**2/(4*n**2)))) / (1 + (z**2/n)) 
#     sd = math.sqrt((p*(1-p)/n)+(z**2/(4*n**2))) / (1 + (z**2/n))
#     p_prime = (p + (z**2/(2*n))) / (1 + (z**2/n))
    return round(score_plus,4), round(score_minus,4)#, round(sd,4), round(p_prime, 4)

parser = argparse.ArgumentParser(description="Plotting script for making the individual graphs "
							"for Supplemental Figure S6.")

parser.add_argument("f", help="Occlusion results csv file, e.g. -f 4eyd_occlusion_totals.csv")
parser.add_argument("prefix", help="Model of glycoform, e.g. 4EYD")

args = parser.parse_args()

gf_names = ["GF 2", "GF 3", "GF 4", "GF 5", "GF 6", "GF 7", "GF 8", "GF 9", "GF 10", "GF 11", "GF 12", "GF 13"]

data_id = []
occlusion_p = []    # list to hold the different proportions for each sample
lag_times = []
independent_samples = []

# -- Open the file and collect the data
inFile = pd.read_csv(args.f)
for i in range(len(inFile.iloc[:,0])):
    data_id.append(inFile.iloc[i,0])
    

# -- collect the proportions
for i in range(len(inFile.iloc[:,4])):
    proportion = round(float(inFile.iloc[i,4]) / 8001,5)
    occlusion_p.append(proportion)
    
for i in range(len(inFile.iloc[:,3])):
    lag_time = float(inFile.iloc[i,3])
    try:
        lag_times.append(lag_time)
        independent_samples.append(round(2000/lag_time,2))
    except:
        lag_times.append(None)
        independent_samples.append(1)

output_file = open(f"output_{args.prefix}_OcclusionBinomial_calculations.csv", 'w')
output_file.write('Sequence_ID, Occlusion_proportion, Independent_Samples, Wilson_Upper, Wilson_Lower\n')

z_alpha = 1.96

wilson_upper = []
wilson_lower = []
for i in range(len(occlusion_p)):
    if occlusion_p[i] == 0:
        wilson_upper.append(0)
        wilson_lower.append(0)
    else:
        up, low = wilsonScore(occlusion_p[i], independent_samples[i], z_alpha)
        # --- This is because the plotting thinks you ADD these on top of where the value 
        # ---    is, instead of where the bounds of the errors are.
        up_adjust = up - occlusion_p[i]
        low_adjust = occlusion_p[i] - low
        # ---
        output_file.write(f"{data_id[i]}, {occlusion_p[i]}, {independent_samples[i]}, {up}, {low}\n")
        wilson_upper.append(up_adjust)
        wilson_lower.append(low_adjust)

occlusion = np.array(occlusion_p)
err_upper = np.array(wilson_upper)
err_lower = np.array(wilson_lower)


output_file.close()

"""
Defining font settings for the figure
"""

rc('font', **{
    'family': 'sans-serif',
    'sans-serif': ['Arial'],
    'size': 10
})
# Set the font used for MathJax - more on this later
rc('mathtext', **{'default': 'regular'})
plt.rc('font', family="Arial")

"""
Code to plot the data
"""
fig, ax = plt.subplots()


#plt.errorbar(gf_names, occlusion, yerr=np.array([err_lower,err_upper]), fmt="ko", markersize=4, ecolor = "red", elinewidth=1.0, capsize=4)
plt.bar(gf_names, occlusion, 0.5, color='lightblue', capsize=3, yerr=np.array([err_lower,err_upper]), ecolor="red")
plt.xticks(rotation=45)
plt.ylim(-0.1, 1.1)
plt.ylabel("Proportion of frames with occlusion")
plt.title(f"Comparison of dimer-interface occlusion proportion: {args.prefix}", weight='bold')
plt.grid(linewidth=0.5)

plt.savefig(f"{args.prefix}_occlusion_binomialPlot.png", dpi=1200)
