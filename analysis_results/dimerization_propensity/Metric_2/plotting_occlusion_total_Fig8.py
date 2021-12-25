#!/usr/bin/env python3

import matplotlib.pyplot as plt 
from matplotlib import rc
import numpy as np


"""
The data used for this figure
"""
acs_names = ["ACS-2", "ACS-3", "ACS-4", "ACS-5", "ACS-6", "ACS-7", "ACS-8", "ACS-9", "ACS-10", "ACS-11", "ACS-12", "ACS-13"]

gf_names = ["GF 2", "GF 3", "GF 4", "GF 5", "GF 6", "GF 7", "GF 8", "GF 9", "GF 10", "GF 11", "GF 12", "GF 13"]

occlude = np.array([7.49906E-05,0,0.116535433,0.995325584,0.603949506,0,0.003274591,0.984901887,0.996150481,0.240819898,0.790526184,0.999475066])
occlude_err = np.array([0.000167684,0,0.038515635,0.004352616,0.192225045,0,0.007322207,0.010475376,0.004494574,0.230783075,0.111019866,0.001173789])



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
fig, ax1 = plt.subplots(figsize=(6,5))

#ax1.errorbar(gf_names, occlude, yerr=occlude_err, fmt="ko", ecolor = "black", capsize=5)
ax1.bar(gf_names, occlude, 0.5, color='lightblue', capsize=3, yerr=occlude_err)
ax1.set_ylim(-0.1, 1.1)
ax1.set_title("Comparison of occlusion proportion across all models", weight='bold')
ax1.grid(linewidth=0.5)

plt.ylabel("Proportion of frames with Occlusion")
plt.xticks(rotation=45)
colors = ['k', 'k', 'k', 'k', 'k', 'k', 'k', 'r', 'r', 'k', 'k', 'r']
for xtick, color in zip(ax1.get_xticklabels(), colors):
    xtick.set_color(color)

plt.savefig("Combined_GFs_comparisonOfOcclusion.png", dpi=600)
