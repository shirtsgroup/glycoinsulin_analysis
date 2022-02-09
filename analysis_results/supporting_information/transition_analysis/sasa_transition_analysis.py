import time
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import rc
from pymbar import timeseries
from prettytable import PrettyTable

def logger(*args, **kwargs):
    print(*args, **kwargs)
    with open('sasa_transition_analysis_results.txt', "a") as f:
        print(file=f, *args, **kwargs)

def analyze_time_series(data):
    [t, g, N_eff] = timeseries.detectEquilibration(data, nskip=1)
    t = 0   # we don't truncate any data here in transition_analysis
    avg = np.mean(data[t:])
    std = np.std(data) / np.sqrt(len(data) / g)

    return avg, std

def analyze_transition(data, f, split):
    data_list = [data[:split], data[split:], data]
    result = []

    for i in data_list:
        avg, std = f(i)
        #print(np.mean(i), np.std(i))
        result.append(f'{avg:.3f} +/- {std:.3f}')

    return result

def tabulate_sasa_result(data_dir, sys, pattern, col, f):
    x = PrettyTable()
    x.add_column('-', ['Change point', 'State A', 'State B', 'Overall'])
    for i in range(len(sys)):
        xvg = data_dir + f'{sys[i]}/' + f'{sys[i].lower()}{pattern}'
        data = np.transpose(np.loadtxt(xvg, comments=["@", "#"]))[col]
        result = [f'{change_loc[i]} ns']
        result.extend(analyze_transition(data, f, change_idx[i]))
        x.add_column(f"{sys[i]}", result)
    
    return x

if __name__ == "__main__":
    t1 = time.time()

    rc("font", **{"family": "sans-serif", "sans-serif": ["DejaVu Sans"], "size": 10})
    # Set the font used for MathJax - more on this later
    rc("mathtext", **{"default": "regular"})
    plt.rc("font", family="Arial")
    
    change_loc = [621.25, 1120.0, 652.5, 1477.5]   # 4EY9, 4EY1, 3I3Z, and 2MVC
    change_idx = [int(change_loc[i] / 0.25) + 1 for i in range(len(change_loc))]  # dt = 250 ps
    sys = ['4EY9', '4EY1', '3I3Z', '2MVC']
    
    logger('[Metric 1: Scissile bond SASA]')
    data_dir = '../../proteolytic_stability/Metric_1/'
    
    logger('(1) SASA of the scissile bond between B25-B26 (units: nm^2)')
    x = tabulate_sasa_result(data_dir, sys, '_sasa_B25.xvg', 2, analyze_time_series)
    logger(x)

    logger('\n(2) SASA of the scissile bond between B25-B26 (units: nm^2)')
    x = tabulate_sasa_result(data_dir, sys, '_sasa_B26.xvg', 2, analyze_time_series)
    logger(x)

    logger('\n[Metric 2: P1 site SASA]')
    data_dir = '../../proteolytic_stability/Metric_2/'
    logger('(1) SASA of residue B24')
    x = tabulate_sasa_result(data_dir, sys, '_sasa_res_B24.xvg', 2, analyze_time_series)
    logger(x)

    logger('\n(2) SASA of residue B25')
    x = tabulate_sasa_result(data_dir, sys, '_sasa_res_B25.xvg', 2, analyze_time_series)
    logger(x)

    t2 = time.time()
    logger(f'Elapsed time: {t2-t1:.2f} ns')    