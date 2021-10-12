#!/usr/bin/env python3

import matplotlib.pyplot as plt 
from matplotlib import rc
import numpy as np
import argparse
import math

def wilsonScore(p, n, z):
    score_plus = (p + (z**2/(2*n)) + z*math.sqrt((p*(1-p)/n) + (z**2/(4*n**2)))) / (1 + (z**2/n))
    score_minus = (p + (z**2/(2*n)) - z*math.sqrt((p*(1-p)/n) + (z**2/(4*n**2)))) / (1 + (z**2/n)) 
#     sd = math.sqrt((p*(1-p)/n)+(z**2/(4*n**2))) / (1 + (z**2/n))
#     p_prime = (p + (z**2/(2*n))) / (1 + (z**2/n))
    return round(score_plus,4), round(score_minus,4)#, round(sd,4), round(p_prime, 4)

parser = argparse.ArgumentParser()

parser.add_argument("f", help="File with total DSSP values for binomial calculation")
parser.add_argument("l", help="File with the Helix lag times for interval calculation")
parser.add_argument("prefix", help="Model of glycoform")

args = parser.parse_args()

gf_names = ["WT", "GF 2", "GF 3", "GF 4", "GF 5", "GF 6", "GF 7", "GF 8", "GF 9", "GF 10", "GF 11", "GF 12", "GF 13"]

data_id = []    # ID for each sample, corresponds to the data in each row
helix_p = []    # list to hold the different proportions for each sample
independent_samples = []

output_file = open(f"output_{args.prefix}_DSSPbinomial_calculations.csv", 'w')
output_file.write('Sequence_ID, Helix_proportion, Independent_Samples, Wilson_Upper, Wilson_Lower\n')

line_count = 1
total = 80010
z_alpha = 1.96
p = 0

with open(args.f, 'r') as f:
    for line in f:
        if line_count == 1:
            line_count += 1
            continue
        else:
            linesplit = line.rstrip("\n").split(",")
            helix_p.append(round(float(linesplit[-1])/total, 5))
            data_id.append(linesplit[0])
            line_count += 1

line_count = 1
with open(args.l, 'r') as f:
    for line in f:
        if line_count == 1:
            line_count += 1
            continue
        else:
            linesplit = line.rstrip("\n").split(",")
            sample = round(2000/float(linesplit[-1]),2)
            independent_samples.append(sample)

wilson_upper = []
wilson_lower = []
for i in range(len(helix_p)):
    up, low = wilsonScore(helix_p[i], independent_samples[i], z_alpha)
    # --- This is because the plotting thinks you ADD these on top of where the value 
    # ---    is, instead of where the bounds of the errors are.
    up_adjust = up - helix_p[i]
    low_adjust = helix_p[i] - low
    # ---
    output_file.write(f"{data_id[i]}, {helix_p[i]}, {independent_samples[i]}, {up}, {low}\n")
    wilson_upper.append(up_adjust)
    wilson_lower.append(low_adjust)

dssp = np.array(helix_p)
dssp_upper = np.array(wilson_upper)
dssp_lower = np.array(wilson_lower)


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

#plt.errorbar(gf_names, dssp, yerr=np.array([dssp_lower,dssp_upper]), fmt="ko", markersize=4, ecolor = "red", elinewidth=1.0, capsize=4)
plt.bar(gf_names, dssp, 0.5, color='lightblue', capsize=3, yerr=np.array([dssp_lower,dssp_upper]), ecolor = "red")
plt.xticks(rotation=45)
plt.ylim(-0.1, 1.1)
plt.ylabel("Proportion of structure with Helix assignment")
plt.title(f"Comparison of Helix Proportion: {args.prefix}", weight='bold')
plt.grid(linewidth=0.5)

plt.savefig(f"{args.prefix}_DSSP_binomialPlot.png", dpi=600)
