import numpy as np 
from prettytable import PrettyTable

def logger(*args, **kwargs):
    print(*args, **kwargs)
    with open('detachment_analysis_results.txt', "a") as f:
        print(file=f, *args, **kwargs)

def get_disorder_fraction(theta_1, theta_2, theta_3):
    # This funciton returns the percentage of each disorder element
    n_b1b7_flip = sum(theta_2 > 10) / len(theta_2) * 100
    n_b1b7_detach = sum(theta_1[theta_2 < 10] > 85) / len(theta_1) * 100
    n_b20b30_detach = sum(theta_3 < 0) / len(theta_3) * 100

    return n_b1b7_flip, n_b1b7_detach, n_b20b30_detach

if __name__ == "__main__":
    models = ['4eyd', '4ey9', '4ey1', '3i3z', '2mvc']
    theta_1, theta_2, theta_3 = [], [], []
    for i in models:
        data = np.transpose(np.loadtxt(f'{i}_detachment.dat'))
        theta_1.append(list(data[1] * 180 / np.pi))  # B1-B7 detachment
        theta_2.append(list(data[2] * 180 / np.pi))  # B1-B7 flipping
        theta_3.append(list(data[3] * 180 / np.pi))  # B20-B30 detachment
    theta_1 = np.array(theta_1).flatten()
    theta_2 = np.array(theta_2).flatten()
    theta_3 = np.array(theta_3).flatten()
    n_frames = len(theta_1)
    s = int(n_frames / len(models))  # size of each trajectory

    x = PrettyTable()
    x.field_names = ['WT model', 'B1-B7 flipping', 'B1-B7 detachment', 'B20-B30 detachment']
    for i in range(5):
        fractions = get_disorder_fraction(theta_1[s * i : s * (i + 1)], theta_2[s * i : s * (i + 1)], theta_3[s * i : s * (i + 1)])
        results = [f'{models[i]}'.upper()]
        results.extend([f'{fractions[i]:.2f}%' for i in range(len(fractions))])
        x.add_row(results)
    
    f = get_disorder_fraction(theta_1, theta_2, theta_3)
    r = ['Total']
    r.extend(f'{f[i]:.2f}%' for i in range(len(f)))
    x.add_row(r)

    logger(x)
