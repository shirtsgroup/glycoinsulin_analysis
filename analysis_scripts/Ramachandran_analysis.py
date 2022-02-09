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

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)


def initialize():
    parser = argparse.ArgumentParser(
        description="This code plot all the Ramachandran plot for the specified residue"
    )
    parser.add_argument("-s", "--sys", help="The system name (in uppercase).")
    parser.add_argument(
        "-r",
        "--res",
        nargs="+",
        type=int,
        help="The residue index of interest (starts from 1).",
    )
    parser.add_argument(
        "-xl",
        "--xlim",
        nargs="+",
        default=[-180, 0],
        help="The minimum and maximum of the phi angle (deg).",
    )
    parser.add_argument(
        "-yl",
        "--ylim",
        nargs="+",
        default=[0, 180],
        help="The minimum and maximum of the psi angle (deg).",
    )
    parser.add_argument(
        "-d",
        "--dir",
        default="./",
        help="The directory of saving the outpus or reading the pickle files, if any.",
    )

    args_parse = parser.parse_args()

    return args_parse


class Logging:
    def __init__(self, output):
        self.output = output

    def logger(self, *args, **kwargs):
        print(*args, **kwargs)
        with open(self.output, "a") as f:
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
    args = initialize()

    t1 = time.time()

    rc("font", **{"family": "sans-serif", "sans-serif": ["DejaVu Sans"], "size": 10})
    # Set the font used for MathJax - more on this later
    rc("mathtext", **{"default": "regular"})
    plt.rc("font", family="serif")

    if args.sys is None:
        path = os.path.abspath(__file__)
        sys_list = ["4EYD", "4EY9", "4EY1", "3I3Z", "2MVC"]
        for i in sys_list:
            if i in path:
                args.sys = i

    if args.sys == "4EYD":
        traj_1 = "../../../wildtype_insulin/4EYD/pH_8.0/Analysis/4eyd_dt250.pdb"
    elif args.sys == "4EY9":
        traj_1 = "../../../wildtype_insulin/4EY9/pH_8.0/Analysis/4ey9_dt250.pdb"
    elif args.sys == "4EY1":
        traj_1 = "../../../wildtype_insulin/4EY1/pH_7.9/Analysis/4ey1_dt250.pdb"
    elif args.sys == "3I3Z":
        traj_1 = "../../../wildtype_insulin/3I3Z/pH_6.9/Analysis/3i3z_dt250.pdb"
    elif args.sys == "2MVC":
        traj_1 = "../../../wildtype_insulin/2MVC/pH_7.3/Analysis/2mvc_dt250.pdb"

    # Step 1: read in save all the trajectories
    # (Note that the minimum version of mdanalysis allow pickling is 2.0.0 beta.)
    L = Logging(f"{args.dir}{args.sys}_Ramachandran_results.txt")
    L.logger(f'\nCommand line: {" ".join(sys.argv)}')
    if os.path.isfile(f"{args.dir}{args.sys.lower()}_trajs.pickle") is True:
        L.logger(f"Reading the trajectories from {args.sys.lower()}_trajs.pickle...")
        with open(f"{args.dir}{args.sys.lower()}_trajs.pickle", "rb") as handle:
            trajs = pickle.load(handle)
    else:
        L.logger(
            f"\nReading in the trajectories and saving {args.sys.lower()}_trajs.pickle ..."
        )
        trajs = []
        u_1 = mda.Universe(traj_1)  # wildtype
        trajs.append(u_1)
        for i in range(2, 14):
            traj_2 = f"glycoform_{i}_ACS/Analysis/glycoform_{i}_ACS_dt250.pdb"
            if os.path.isfile(traj_2):
                u_2 = mda.Universe(traj_2)  # glycoform
            else:
                u_2 = "missing"
            trajs.append(u_2)
        with open(f"{args.dir}{args.sys.lower()}_trajs.pickle", "wb") as handle:
            pickle.dump(trajs, handle, protocol=pickle.HIGHEST_PROTOCOL)

    for r in args.res:
        sec_str = f"\nStart analyzing residue {r} ..."
        L.logger(sec_str)
        L.logger("=" * len(sec_str))
        suffix = f"res_{r}"  # for the output filename
        title_str = f"residue B{r - 21}"

        if os.path.isfile(f"{args.dir}{args.sys.lower()}_rama_res{r}.pickle") is True:
            L.logger(
                f"Reading the Ramachandran plot data from {args.sys.lower()}_rama_res{r}.pickle ..."
            )
            with open(
                f"{args.dir}{args.sys.lower()}_rama_res{r}.pickle", "rb"
            ) as handle:
                rama_data = pickle.load(handle)
        else:
            L.logger(
                f"Generating Ramchandran plot data and saving as {args.sys.lower()}_rama_res{r}.pickle ..."
            )
            rama_data = []

            u_1 = trajs[0]
            region_1 = u_1.select_atoms(f"resid {r}")
            L.logger(f'Residue of interest: {region_1.resnames[0]} B{r - 21}') # Just print once
            rama_1 = dihedrals.Ramachandran(region_1).run()
            rama_data.append(rama_1)

            for i in range(2, 14):
                if trajs[i - 1] == "missing":
                    rama_data.append("missing")
                else:
                    u_2 = trajs[i - 1]
                    region_2 = u_2.select_atoms(f"resid {r}")
                    rama_2 = dihedrals.Ramachandran(region_2).run()
                    rama_data.append(rama_2)
            with open(
                f"{args.dir}{args.sys.lower()}_rama_res{r}.pickle", "wb"
            ) as handle:
                pickle.dump(rama_data, handle, protocol=pickle.HIGHEST_PROTOCOL)

        # Step 2: Plot the Ramachandran plots!
        L.logger("Plotting Ramachandran plots ...")
        plt.figure(figsize=(12, 8))
        for i in range(2, 14):
            if trajs[i - 1] == "missing":
                plt.subplot(3, 4, i - 1)
                plt.text(
                    -175,
                    5,
                    f"* GF {i} data missing",
                    fontsize=7,
                    weight="bold",
                )
                rama_data[0].plot(
                    color="yellow", marker=".", ref=True, alpha=0.5, label="WT", s=0.1
                )
            else:
                plt.subplot(3, 4, i - 1)
                rama_data[0].plot(
                    color="yellow", marker=".", ref=True, alpha=0.5, label="WT", s=0.1
                )
                rama_data[i - 1].plot(
                    color="red", marker=".", alpha=0.5, label=f"GF {i}", s=0.1
                )

            plt.arrow(-130, 110, -10, -15, head_width=3, head_length=3, ec='black', fc='black')
            plt.text(-160, 68, "β-sheet\nregion", fontsize=7)

            if args.xlim is not None:
                plt.xlim(args.xlim[0], args.xlim[1])
            if args.ylim is not None:
                plt.ylim(args.ylim[0], args.ylim[1])
            plt.grid()

            lgnd = plt.legend(loc="lower right", scatterpoints=1, fontsize=10)
            for handle in lgnd.legendHandles:
                handle.set_sizes([30.0])
            plt.xlabel("$\phi$ (deg)")
            plt.ylabel("$\psi$ (deg)")

        plt.suptitle(
            f"\nThe Ramachandran plots of {title_str} of all {args.sys}-based glycoforms",
            weight="bold",
            fontsize=15,
        )
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig(f"{args.dir}{args.sys}_multi_rama_plot_{suffix}.png", dpi=600)

        # Step 3: Calculate and plot the percentage of the points in the beta-sheet region
        L.logger(
            f"Calculating the percentage of the poitns in the beta-sheet region ...\n"
        )
        RA = RamachandranAnalyzer()
        beta_fractions, sys, missing = [], [], []
        for i in range(len(rama_data)):
            if i == 0:
                sys.append("WT")
            else:
                sys.append(f"GF {i + 1}")

            if rama_data[i] == "missing":
                beta_fractions.append(0)
                missing.append(sys[i])
                L.logger(f"The trajectory data of {sys[i]} is missing.")
            else:
                angles = rama_data[i].angles
                n1, n2, n3 = angles.shape  # n1: n_frames, n2: n_residues, n3=2
                angles = angles.reshape(n1 * n2, n3)
                f = RA.beta_sheet_fraction(angles)
                beta_fractions.append(f)
                L.logger(
                    f"Beta-sheet fraction of residue B{r - 21} in {sys[i]}: {f:.2f}%"
                )
        with open(f"{args.dir}{args.sys.lower()}_beta_res{r}.pickle", "wb") as handle:
            pickle.dump(beta_fractions, handle, protocol=pickle.HIGHEST_PROTOCOL)

        exp = ['$=$', '$-$', '$-$', '$-$', r'$\approx$', '$-$', '$-$', '$-$', r'$\approx$', '$+$', '$-$', r'$\approx$', '$+$']
        fig, ax = plt.subplots(figsize=(12, 3))
        plt.bar([f'{sys[i]}\n{exp[i]}' for i in range(len(sys))], beta_fractions, color="lightblue")
        plt.xticks(fontsize=7)
        plt.ylabel("Percentage of the \nframes in the β-sheet region (%)")
        if len(missing) != 0:
            plt.text(
                0.02,
                0.92,
                f"*Missing data: {', '.join(missing)}",
                transform=ax.transAxes,
                fontsize=7,
                weight="bold",
            )
        plt.title(f"{args.sys} (Residue B{r - 21})", weight="bold")
        plt.grid()
        plt.savefig(f"{args.dir}{args.sys}_beta_fraction_{suffix}.png", dpi=600)

    t2 = time.time()
    L.logger(f"\nElapsed time: {t2 - t1} seconds.")

