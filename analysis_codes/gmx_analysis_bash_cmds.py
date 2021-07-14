import os
import argparse
import glob
import gromacs as gmx
import MDAnalysis


def initialize():
    parser = argparse.ArgumentParser(
        description="This code use GROMACS wrapper to carry out GROMACS commands for data analysis."
    )
    parser.add_argument(
        "-s",
        "--sys",
        required=True,
        help="System name, the prefix of the input files.",
    )
    parser.add_argument("-n", "--ndx", help="The filename of the NDX input.")
    parser.add_argument("-g", "--gro", help="The filename of the GRO input.")
    parser.add_argument("-t", "--tpr", help="The filename of the TPR input.")
    parser.add_argument("-x", "--xtc", help="The filename of the XTC input.")
    parser.add_argument(
        "-sk", "--skip", type=int, default=0, help="The number of steps to be skipped."
    )

    args_parse = parser.parse_args()
    return args_parse


def get_gmx_fn(module, func_name):
    """
    func_name (str): The name of the method in GROMACS wrapper
    """
    try:
        return getattr(module, func_name)
    except AttributeError:
        func_name += "_mpi"
        return getattr(module, func_name)


if __name__ == "__main__":
    args = initialize()

    # Step 1: Setting things up
    if args.gro is None:
        args.gro = f"../MD/{args.sys}_md.gro"
    if args.ndx is None:
        args.ndx = f"{args.sys}.ndx"
    if args.tpr is None:
        args.tpr = f"../MD/{args.sys}_md.tpr"
    if args.xtc is None:
        args.xtc = f"../MD/{args.sys}_md.xtc"

    # Step 2: Create the following index groups:
    if not glob.glob("*.ndx"):  # no ndx file
        os.system(f'mpirun -np 1 gmx_mpi make_ndx -f {args.gro} -o {args.ndx}')

    ndx = gmx.fileformats.ndx.NDX()
    ndx.read(args.ndx)
    n_groups = len(ndx)

    # Add index group (1): CONH atoms between B25 (index 45) and B26 (index 46)
    u = MDAnalysis.Universe(args.gro)
    C_idx = u.residues[45].atoms[-2].index + 1
    ndx["B25-B26"] = [C_idx, C_idx + 1, C_idx + 2, C_idx + 3]

    # Add index group (2): CONH atoms between B26 (index 45) and B27 (index 46)

    C_idx = u.residues[46].atoms[-2].index + 1
    ndx["B26-B27"] = [C_idx, C_idx + 1, C_idx + 2, C_idx + 3]

    # Add index group (3): The whole glycoform
    for i in range(len(u.residues) - 1):
        if u.residues[i].resname != "SOL" and u.residues[i + 1].resname == "SOL":
            res_idx = i
            break

    glycoform_atoms_idx = list(u.residues[: res_idx + 1].atoms.indices + 1)
    ndx["glycoform"] = glycoform_atoms_idx
    ndx.write(args.ndx)

    # Step 3: Data analysis
    # Calcualte SASA calculation (note that the index in the .ndx file also starts from 0)
    groups = list(ndx.keys())
    non_water_idx = groups.index("non-Water")
    B25_idx = groups.index("B25-B26")
    B26_idx = groups.index("B26-B27")

    if args.skip >= 1:
        pass
    else:
        os.system(f'mpirun -np 1 gmx_mpi sasa -f {args.xtc} -s {args.tpr} -n {args.ndx} -o {args.sys}_sasa_B25.xvg -surface "group {non_water_idx}" -output "group {B25_idx}" -dt 250')

    if args.skip >= 2:
        pass
    else:
        os.system(f'mpirun -np 1 gmx_mpi sasa -f {args.xtc} -s {args.tpr} -n {args.ndx} -o {args.sys}_sasa_B25.xvg -surface "group {non_water_idx}" -output "group {B26_idx}" -dt 250')

    # Convert the trajectory for H.B. analysis
    os.system(f'echo 1 | mpirun -np 1 gmx_mpi trjconv -f {args.xtc} -s {args.tpr} -n {args.ndx} -o {args.sys}_dt250.pdb -dt 250')
