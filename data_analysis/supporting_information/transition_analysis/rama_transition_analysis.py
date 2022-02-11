import time
import numpy as np
import warnings
import MDAnalysis as mda
import matplotlib.pyplot as plt
from matplotlib.path import Path
from MDAnalysis.analysis import dihedrals
from MDAnalysis.analysis.data.filenames import Rama_ref
from prettytable import PrettyTable

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

def logger(*args, **kwargs):
    print(*args, **kwargs)
    with open('rama_transition_analysis_results.txt', "a") as f:
        print(file=f, *args, **kwargs)

class RamachandranAnalyzer:
    def __init__(self):
        """
        The init functions of the class defines the coordinates enclosing the beta-sheet region
        in a Ramachandran plot. The reference was taken from the source code of Ramachandran.plot
        from MDAnalysis (see https://shorturl.at/fnBCN).

        Attributes
        ----------
        beta_region (matplotlib.path.Path): matplotlib Path object defining the region of
                                            beta-sheets on a Ramachandran plot.
        """
        X, Y = np.meshgrid(np.arange(-180, 180, 4), np.arange(-180, 180, 4))
        # levels are set below so the contour line corresponding to the β-sheet region will be shown
        levels = [1, 17, 15000]
        cs = plt.contourf(X, Y, np.load(Rama_ref), levels=levels)
        coords = cs.collections[0].get_paths()[34].vertices  # coordaintes
        # x and y below roughly defines the beta-sheet region
        x = np.transpose(coords)[0]
        y = np.transpose(coords)[1]

        # Slice the data to better sketch the region
        x = x[(y > 60)][50:160]
        y = y[(y > 60)][50:160]
        additional_x = np.arange(x[-1], x[0])
        additional_y = np.ones(len(additional_x)) * 180
        x = np.concatenate((x, additional_x))
        y = np.concatenate((y, additional_y))

        # Turn x, y into a list of tuples and make it as a path
        verts = []
        for i in range(len(x)):
            verts.append((x[i], y[i]))
        self.beta_region = Path(verts)

    def beta_sheet_fraction(self, rama_points):
        """
        Parameters
        ----------
        rama_points (np.ndarray): 
            An N by 2 numpy array containing the coordinates of the points to be checked.

        Returns
        -------
        beta_frac (float): The percentage of the points in the beta-sheet region.
        """
        grid = self.beta_region.contains_points(rama_points)  # len: n_frrames
        beta_fractions = np.sum(grid) / len(grid) * 100

        return beta_fractions


if __name__ == "__main__":
    t1 = time.time()

    trajs = []  # working directory: /ocean/projects/cts160011p/wehs7661/Glycoinsulin_project/wildtype_insulin
    trajs.append("4EY9/pH_8.0/Analysis/4ey9_dt250.pdb")
    trajs.append("4EY1/pH_7.9/Analysis/4ey1_dt250.pdb")
    trajs.append("3I3Z/pH_6.9/Analysis/3i3z_dt250.pdb")
    trajs.append("2MVC/pH_7.3/Analysis/2mvc_dt250.pdb")

    sys = ['4EY9', '4EY1', '3I3Z', '2MVC']

    change_loc = [621.25, 1120.0, 652.5, 1477.5]   # 4EY9, 4EY1, 3I3Z, and 2MVC
    change_idx = [int(change_loc[i] / 0.25) + 1 for i in range(len(change_loc))]  # dt = 250 ps

    logger('[Beta-sheet propensity for each state of 4 wild-type insulin models]')
    RA = RamachandranAnalyzer()
    for i in range(len(trajs)):
        logger(f'({i + 1}) Wild-type model: {sys[i]} (Change point: {change_loc[i]} ns)')
        # One table for each WT variant
        x = PrettyTable()
        x.add_column('-', ['State A', 'State B', 'Overall'])

        u = mda.Universe(trajs[i])
        for j in range(43, 47):   # residues 43 to 47
            region = u.select_atoms(f"resid {j}")
            rama = dihedrals.Ramachandran(region).run()
            angles = rama.angles
            n1, n2, n3 = angles.shape  # n1: n_frames, n2: n_residues, n3=2
            angles = angles.reshape(n1 * n2, n3)
            angles_A = angles[:change_idx[i]]
            angles_B = angles[change_idx[i]:]

            # Calculate beta-sheet propensity
            f = RA.beta_sheet_fraction(angles)  
            f_A = RA.beta_sheet_fraction(angles_A)
            f_B = RA.beta_sheet_fraction(angles_B)

            x.add_column(f'Residue B{j - 21}', [f'{f_A:.2f}%', f'{f_B:.2f}%', f'{f:.2f}%'])
        
        logger(x)
        logger('\n')

    t2 = time.time()
    logger(f'Time elapsed: {t2 - t1:.2f} seconds.')
