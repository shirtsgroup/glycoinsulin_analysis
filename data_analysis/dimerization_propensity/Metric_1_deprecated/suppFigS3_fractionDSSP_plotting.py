import matplotlib.pyplot as plt 
from matplotlib import rc
import numpy as np
import argparse

parser = argparse.ArgumentParser(description="Plotting script to plot the individual graphs"
								" used to generate Supplemental Figure S3.")

parser.add_argument("f", help="File with total DSSP values for fraction structure plotting.")
parser.add_argument("p", help="Model of glycoform")

args = parser.parse_args()

gf_names = ["WT", "GF 2", "GF 3", "GF 4", "GF 5", "GF 6", "GF 7", "GF 8", "GF 9", "GF 10", "GF 11", "GF 12", "GF 13"]

coil = []
bend = []
turn = []
helix = []

line_count = 1
total = 80010
with open(args.f, 'r') as f:
    for line in f:
        if line_count == 1:
            line_count += 1
            continue
        else:
            linesplit = line.rstrip("\n").split(",")
            coil.append(round(int(linesplit[2])/total, 5))
            bend.append(round(int(linesplit[3])/total, 5))
            turn.append(round(int(linesplit[4])/total, 5))
            helix.append(round(int(linesplit[5])/total, 5))
            line_count += 1

coil = np.array(coil)
bend = np.array(bend)
turn = np.array(turn)
helix = np.array(helix)

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
fig, ax = plt.subplots(figsize=(7,5))
width = 0.45

plt.bar(gf_names, coil, width, fill=True, facecolor = 'White', edgecolor="black", hatch="/", label="Coil")
plt.bar(gf_names, bend, width, bottom=coil, fill=True, facecolor="dimgrey", edgecolor="black", label="Bend")
plt.bar(gf_names, turn, width, bottom=(coil+bend), fill=True, facecolor='white', edgecolor="black", hatch="++", label="Turn")
plt.bar(gf_names, helix, width, bottom=(coil+bend+turn), fill=True, facecolor="lightgrey", edgecolor="black", label="Helix")

ax.set_axisbelow(True)
plt.xticks(rotation=45)
plt.ylim(0, 1.1)
plt.ylabel("Fraction secondary structure")
plt.title(f"DSSP secondary structure analysis: {args.p} LeuB17-TyrB26", weight='bold')
plt.legend(framealpha=1, edgecolor="black", fontsize="small", bbox_to_anchor=(0.96, 0.22), loc="upper left")
plt.grid(linewidth=0.5)

# plt.show()
plt.savefig(f"{args.p}_dssp.png", dpi=1200)
