import numpy as np
import matplotlib.pyplot as plt 
from matplotlib import rc 
from colour import Color

if __name__ == '__main__':
    exp = ['$-$', '$-$', '$-$', '$≈$', '$-$', '$-$', '$-$', '$≈$', '$+$', '$-$', '$≈$', '$+$']
    sys = [f"GF {i}" for i in range(2, 14)]
    labels = [f'{sys[i]}\n({exp[i]})' for i in range(len(sys))]

    txt_colors = []
    for i in range(len(labels)):
        if "$+$" in labels[i]:
            txt_colors.append('blue')
        elif r"$-$" in labels[i]:
            txt_colors.append('red')
        else:
            txt_colors.append('magenta')

    height = []
    height.append([0])              # GF 2
    height.append([0])              # GF 3
    height.append([0])              # GF 4
    height.append([16.8])           # GF 5
    height.append([12.8])           # GF 6
    height.append([0])              # GF 7
    height.append([0])              # GF 8
    height.append([9.8, 6.4])       # GF 9
    height.append([41.6, 32.4, 8])  # GF 10
    height.append([0])              # GF 11
    height.append([32.2])           # GF 12
    height.append([41.2, 25, 9.8, 5.8])   # GF 13

    err = []
    err.append([0])                 # GF 2
    err.append([0])                 # GF 3
    err.append([0])                 # GF 4
    err.append([16.36])             # GF 5
    err.append([11.63])             # GF 6
    err.append([0])                 # GF 7
    err.append([0])                 # GF 8
    err.append([12.11, 8.06])       # GF 9
    err.append([9.33, 12.11, 6.66]) # GF 10
    err.append([0])                 # GF 11
    err.append([10.57])             # GF 12
    err.append([11.91, 12.79, 5.04, 7.11])  # GF 13

    pos = []
    pos.append([0])                 # GF 2
    pos.append([1])                 # GF 3
    pos.append([2])                 # GF 4
    pos.append([3])                 # GF 5
    pos.append([4])                 # GF 6
    pos.append([5])                 # GF 7
    pos.append([6])                 # GF 8
    pos.append([6.875, 7.125])      # GF 9
    pos.append([7.775, 8, 8.225])   # GF 10
    pos.append([9])                 # GF 11
    pos.append([10])                # GF 12
    pos.append([10.625, 10.875, 11.125, 11.375])

    hbond = []
    hbond.append(['No H-bond'])         # GF 2
    hbond.append(['No H-bond'])         # GF 3
    hbond.append(['No H-bond'])         # GF 4
    hbond.append(['Phe(B24)-GalNAc(1)'])  # GF 5
    hbond.append(['Thr(B27)-GalNAc(1)'])  # GF 6
    hbond.append(['No H-bond'])         # GF 7
    hbond.append(['No H-bond'])         # GF 8
    hbond.append(['Phe(B24)-Man(1)', 'Tyr(B16)-Man(1)'])    # GF 9
    hbond.append(['Thr(B27)-Man(2)', 'Phe(B24)-Man(1)', 'Tyr(B16)-Man(1)'])   # GF 10
    hbond.append(['No H-bond'])         # GF 11
    hbond.append(['Thr(B30)-Man(2)'])   # GF 12
    hbond.append(['Thr(B27)-Man(2)', 'Phe(B24)-Man(1)', 'Tyr(B16)-Man(1)', 'Thr(B27)-Man(3)'])

    all_heights, all_err, all_pos, all_hbonds = [], [], [], []
    for i in range(len(height)):
        for j in range(len(height[i])):
            all_heights.append(height[i][j])
            all_err.append(err[i][j])
            all_pos.append(pos[i][j])
            all_hbonds.append(hbond[i][j])

    # assign colors to different types of hbonds
    colors = []
    for i in range(len(all_hbonds)):
        if 'Phe(B24)' in all_hbonds[i] or 'Thr(B27)' in all_hbonds[i]:
            colors.append('blue')  # could be other color as needed 
        else:
            colors.append('black')

    rc('font', **{
       'family': 'sans-serif',
       'sans-serif': ['DejaVu Sans'],
       'size': 10
    })
    # Set the font used for MathJax - more on this later
    rc('mathtext', **{'default': 'regular'})
    plt.rc('font', family='Arial')

    fig = plt.figure(figsize=(15, 6))
    ax1 = fig.add_subplot(111)
    """
    ax2 = ax1.twiny()
    fig.subplots_adjust(bottom=0.2)   # extra space for the 2nd axis
    """

    for i in range(len(all_heights)):
        #if all_heights[i] == 0:
        #    ax1.text(all_pos[i] - 0.18, 1.5, 'None', fontsize=12)
        ax1.bar(all_pos[i], all_heights[i], 0.2, color='lightblue', capsize=3, yerr=all_err[i])

    pos_1 = np.arange(12) - 0.2
    pos_2 = np.arange(12) - 0.1
    for i in range(len(sys)):
        ax1.text(pos_1[i], 52, f'{sys[i]}', fontsize=15, weight='bold', c=txt_colors[i])
        ax1.text(pos_2[i], 48, f'({exp[i]})', fontsize=15, weight='bold', c=txt_colors[i])

    # plot the grid lines
    h_pos = np.linspace(-5, 55, 26)
    v_pos = np.linspace(0, 11, 12) + 0.5
    for i in range(len(h_pos)):
        ax1.axhline(y=h_pos[i], color='lightgrey', linewidth=0.5)
    for i in range(len(v_pos)):
        ax1.axvline(x=v_pos[i], color='black', linestyle='--', alpha=0.7, linewidth=0.5)
    
    ax1.set_xticks(all_pos)
    ax1.set_xticklabels(all_hbonds, rotation=50, fontsize=16, ha='right')
    for color, tick in zip(colors, ax1.xaxis.get_major_ticks()):
        tick.label1.set_color(color)
    ax1.set_xlim([-0.5, 11.5])
    ax1.set_ylabel('Average existence percentage (%)', fontsize=16)
    """
    ax2.set_xlim([-0.5, 10.5])
    new_ticks = np.linspace(0, 11, 12) 
    ax2.xaxis.set_ticks_position("bottom")
    ax2.xaxis.set_label_position("bottom")
    ax2.set_frame_on(True)
    ax2.patch.set_visible(False)
    ax2.set_xticks(new_ticks)
    ax2.spines["bottom"].set_position(("axes", -0.25))
    ax2.set_xticklabels(labels, fontsize=18, weight='bold')
    """
    plt.title('Glycan-involved hydrogen bonds of all glycoforms', weight='bold', fontsize=16)
    plt.tight_layout()
    plt.savefig('hbond_results', dpi=600)