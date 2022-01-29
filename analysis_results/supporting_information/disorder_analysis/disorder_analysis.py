import numpy as np 

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
    n_frames = len(theta_1)  # same for all time series

    # B1-B7 flipping 
    no_flip_idx = []
    for i in range(len(theta_2)):
        if theta_2[i] < 10:
            no_flip_idx.append(i)
        
    # B1-B7 detachment
    n_b1b7_detach = 0
    for i in no_flip_idx:
        if theta_1[i] > 85:
            n_b1b7_detach += 1

    # B20-B30 detachment
    n_b20b30_detach = 0
    for i in theta_3:
        if i < 0:
            n_b20b30_detach += 1

    print(f'{n_b1b7_detach / n_frames * 100:.2f}% of the insulin structures have B1-B7 detachment.')    
    print(f'{n_b20b30_detach / n_frames * 100:.2f}% of the insulin structures have B20-B30 detachment.')    


    