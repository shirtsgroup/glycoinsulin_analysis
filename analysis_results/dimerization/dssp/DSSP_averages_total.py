#!/usr/bin/env python3

import matplotlib.pyplot as plt 
from matplotlib import rc
import numpy as np


"""
The data used for this figure
"""
acs_names = ["ACS-2", "ACS-3", "ACS-4", "ACS-5", "ACS-6", "ACS-7", "ACS-8", "ACS-9", "ACS-10", "ACS-11", "ACS-12", "ACS-13"]

gf_names = ["GF 2", "GF 3", "GF 4", "GF 5", "GF 6", "GF 7", "GF 8", "GF 9", "GF 10", "GF 11", "GF 12", "GF 13"]

helix = np.array([0.1180,0.1409,0.0833,0.0765,0.1594,0.0952,0.0559,0.1382,0.0981,0.0968,0.0524,0.1372])
helix_err = np.array([0.0425,0.0406,0.0489,0.0479,0.0533,0.0803,0.0479,0.0413,0.0649,0.0685,0.0498,0.0610])

turn = np.array([0.1548,0.1313,0.1711,0.1816,0.1427,0.1566,0.1942,0.1426,0.1545,0.1549,0.1805,0.1344])
turn_err = np.array([0.0258,0.0218,0.0240,0.0291,0.0226,0.0389,0.0478,0.0275,0.0419,0.0399,0.0187,0.0329])

bend = np.array([0.1267,0.1291,0.1421,0.1458,0.1172,0.1530,0.1492,0.1207,0.1406,0.1444,0.1647,0.1266])
bend_err = np.array([0.0195,0.0176,0.0255,0.0384,0.0219,0.0315,0.0357,0.0171,0.0265,0.0283,0.0245,0.0222])

coil = np.array([0.6005,0.5988,0.6036,0.5961,0.5807,0.5951,0.6006,0.5985,0.6068,0.6040,0.6024,0.6017])
coil_err = np.array([0.0042,0.0036,0.0014,0.0232,0.0403,0.0161,0.0030,0.0034,0.0059,0.0042,0.0082,0.0061])


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
fig, ax1 = plt.subplots()

#ax1.errorbar(gf_names, helix, yerr=helix_err, fmt="ko", ecolor = "black", capsize=5)
ax1.bar(gf_names, helix, 0.5, color='lightblue', capsize=3, yerr=helix_err, ecolor="black")
ax1.set_ylim(0.00, 0.25)
ax1.set_title("Comparison of Helix Assignments", weight='bold')
ax1.set_ylabel("Fraction of Secondary Structure")
ax1.grid(linewidth=0.5)
plt.xticks(rotation=45)

plt.savefig("Combined_GFs_helixAssignmentOnly.png", dpi=600)


fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=True, figsize=(13,9))

ax1.bar(gf_names, turn, 0.5, color='lightblue', capsize=3, yerr=turn_err, ecolor="black")
ax1.set_ylim(0.00, 0.25)
ax1.set_title("Comparison of Turn Assignments", weight='bold')
ax1.grid(linewidth=0.5)

ax2.bar(gf_names, bend, 0.5, color='lightblue', capsize=3, yerr=bend_err, ecolor="black")
ax2.set_ylim(0.00, 0.25)
ax2.set_title("Comparison of Bend Assignments", weight='bold')
ax2.grid(linewidth=0.5)

ax3.bar(gf_names, coil, 0.5, color='lightblue', capsize=3, yerr=coil_err, ecolor="black")
ax3.set_ylim(0.45, 0.70)
ax3.set_title("Comparison of Coil Assignments", weight='bold')
ax3.grid(linewidth=0.5)

ax4.set_axis_off()

fig.supylabel("Fraction of Secondary Structure", x=0.065, fontsize='x-large')
plt.xticks(rotation=45)
plt.sca(ax3)
plt.xticks(rotation=45)

plt.savefig("Combined_GFs_other3DSSP.png", dpi=600)
