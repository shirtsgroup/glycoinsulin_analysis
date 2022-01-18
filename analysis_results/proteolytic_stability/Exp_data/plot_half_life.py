import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc


def read_experimental_data():
    f = open("experimental_data.txt")
    lines = f.readlines()
    f.close()

    data, height, err = [], [], []
    for i in lines:
        data.append(float(i.split(",")[-1]))

    for i in range(len(data)):
        if i % 2 == 0:  # height
            height.append(data[i])
        else:
            err.append(data[i] - data[i - 1])

    return height, err


if __name__ == "__main__":
    rc("font", **{"family": "sans-serif", "sans-serif": ["DejaVu Sans"], "size": 10})
    # Set the font used for MathJax - more on this later
    rc("mathtext", **{"default": "regular"})
    plt.rc("font", family="serif")

    sys = []
    sys.append("Wild type")
    sys.append(r"GalNAc$\alpha$-SerA9")
    sys.append(r"GalNAc$\alpha$-SerA12")
    sys.append(r"GalNAc$\alpha$-SerB9")
    sys.append(r"GalNAc$\alpha$-ThrB27")
    sys.append(r"GalNAc$\alpha$-ThrB30")
    sys.append(r"Man$\alpha$-SerA9")
    sys.append(r"Man$\alpha$2Man$\alpha$-SerA9")
    sys.append(r"Man$\alpha$-ThrB27")
    sys.append(r"Man$\alpha$2Man$\alpha$-ThrB27")
    sys.append(r"Man$\alpha$-ThrB30")
    sys.append(r"Man$\alpha$2Man$\alpha$-ThrB30")
    sys.append(r"Man$\alpha$2Man$\alpha$2Man$\alpha$-ThrB27")

    h, err = read_experimental_data()

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()
    fig.subplots_adjust(bottom=0.2)   # extra space for the 2nd axis

    color = 'white'

    ax1.bar(
        range(1, 14),
        h,
        width=0.8,
        yerr=err,
        error_kw=dict(ecolor="red", capsize=3),
        color="lightblue",
    )
    ax1.spines["bottom"].set_position(("axes", -0.12))
    ax1.set_xlim([0, 14])
    ax1.set_xticks(range(1, 14))
    ax1.tick_params(axis="x", colors=color)
    ax1.tick_params(axis="y", colors=color)
    ax1.set_xticklabels(sys, rotation=50, ha="right", color=color)
    ax1.set_ylabel(r"$\alpha$-Chymotrypsin half-life (min)", color=color)
    ax1.tick_params(axis=u'both', which=u'both',length=0)

    for spine in ax1.spines.values():
        spine.set_edgecolor(color)
    
    ax2.set_xticklabels([str(i) for i in range(1, 14)], weight='bold', color=color)
    ax2.set_frame_on(True)
    ax2.patch.set_visible(False)
    ax2.set_xlim([0, 14])
    ax2.set_xticks(range(1, 14))
    ax2.tick_params(axis="x", colors=color)
    ax2.xaxis.set_ticks_position("bottom")
    ax2.xaxis.set_label_position("bottom")

    for spine in ax2.spines.values():
        spine.set_edgecolor(color)
    
    ax1.grid()
    ax1.set_title('Experimental work by Guan et al.', weight='bold', color=color)
    plt.tight_layout()
    plt.savefig("half_life.png", dpi=600, transparent=True)
