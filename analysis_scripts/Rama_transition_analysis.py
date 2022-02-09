import argparse
import os
import pickle
import sys
import time
import warnings

import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
from matplotlib import rc
from matplotlib.path import Path
from MDAnalysis.analysis import dihedrals
from MDAnalysis.analysis.data.filenames import Rama_ref

def logger(*args, **kwargs):
    print(*args, **kwargs)
    with open('Rama_transition_analysis_results.txt', "a") as f:
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

    def beta_sheet_fraction(self, rama_obj):
        """
        Parameters
        ----------
        rama_obj (MDAnalysis object):
            An MDAnalysis.analysis.dihedrals.Ramachandran object generated 
            by dihedrals.Ramachandran(region).run()

        Returns
        -------
        beta_frac (float): The percentage of the points in the beta-sheet region.
        """
        angles = rama_obj.angles
        n1, n2, n3 = angles.shape
        angles = angles.reshape(n1* n2, n3)

        grid = self.beta_region.contains_points(angles)  # len: n_frrames
        beta_fractions = np.sum(grid) / len(grid) * 100

        return beta_fractions



if __name__ == "__main__":
    t1 = time.time()

    trajs = []
    trajs.append("../../../wildtype_insulin/4EY9/pH_8.0/Analysis/4ey9_dt250.pdb")
    trajs.append("../../../wildtype_insulin/4EY1/pH_7.9/Analysis/4ey1_dt250.pdb")
    trajs.append("../../../wildtype_insulin/3I3Z/pH_6.9/Analysis/3i3z_dt250.pdb")
    trajs.append("../../../wildtype_insulin/2MVC/pH_7.3/Analysis/2mvc_dt250.pdb")

    state_A, state_B, overall = [], [], []
    RA = RamachandranAnalyzer()
    for i in range(len(trajs)):
        u = mda.Universe(trajs[i])
        for j in range(43, 47):   # residues 43 to 47
            region = u.select_atoms(f"resid {j}")
            rama = dihedrals.Ramachandran(region).run()
            angles = rama.angles
            f = RA.beta_sheet_fraction(rama)  # beta-sheet propensity

    t2 = time.time()
