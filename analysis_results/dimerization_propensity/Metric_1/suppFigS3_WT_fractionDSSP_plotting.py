import matplotlib.pyplot as plt 
from matplotlib import rc
import numpy as np

"""

Plotting for: DAR1-10a 4eyd DSSP plots

"""

"""
The data used for this figure
"""

gf_names = ["4EYD", "4EY1", "4EY9", "3I3Z", "2MVC"]

coil = np.array([0.60, 0.60, 0.60, 0.60, 0.60])
bend = np.array([0.13, 0.14, 0.10, 0.10, 0.15])
turn = np.array([0.15, 0.15, 0.13, 0.09, 0.20])
helix = np.array([0.12, 0.11, 0.17, 0.21, 0.05])


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

fig, ax = plt.subplots(figsize=(6.5,5))
width = 0.45
ax.set_axisbelow(True)

plt.bar(gf_names, coil, width, fill=True, facecolor = 'white', edgecolor="black", hatch='/', label="Coil")
plt.bar(gf_names, bend, width, bottom=coil, fill=True, facecolor="dimgrey", edgecolor="black", label="Bend")
plt.bar(gf_names, turn, width, bottom=(coil+bend), fill=True, facecolor='white', edgecolor="black", hatch="++", label="Turn")
plt.bar(gf_names, helix, width, bottom=(coil+bend+turn), fill=True, facecolor="lightgrey", edgecolor="black", label="Helix")

plt.xticks(rotation=45)
plt.ylim(0, 1.1)
plt.ylabel("Fraction secondary structure")
plt.title("DSSP secondary structure analysis: wild type models", weight='bold')
plt.legend(framealpha=1, edgecolor="black", fontsize="small", bbox_to_anchor=(0.96, 0.22), loc="upper left")
plt.grid(linewidth=0.5)

# plt.show()
plt.savefig("WT_overall_dssp.png", dpi=600)
