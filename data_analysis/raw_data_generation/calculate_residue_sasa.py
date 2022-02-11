import os
import argparse
import glob
import gromacs as gmx
import MDAnalysis


def initialize():
    parser = argparse.ArgumentParser(
        description="This code uses GROMACS wrapper to carry out GROMACS commands for data analysis."
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
        args.xtc = glob.glob("../MD/*.xtc")[0]

    # Step 2: Create the following index groups:
    if not glob.glob("*.ndx"):  # no ndx file
        os.system(f'mpirun -np 1 gmx_mpi make_ndx -f {args.gro} -o {args.ndx}')

    ndx = gmx.fileformats.ndx.NDX()
    ndx.read(args.ndx)
    n_groups = len(ndx)

    # Add index group (1): residue B24 (index 44, the whole residue)
    u = MDAnalysis.Universe(args.gro)
    if "res-B24" not in list(ndx.keys()):
        ndx["res-B24"] = list(u.residues[44].atoms.indices + 1)

    # Add index group (2): residue B25 (index 45, the whole residue)
    if "res-B25" not in list(ndx.keys()):
        ndx["res-B25"] = list(u.residues[45].atoms.indices + 1)

    ndx.write(args.ndx)

    # Step 3: SASA calculation
    # Note that the index in the .ndx file also starts from 0)
    groups = list(ndx.keys())
    non_water_idx = groups.index("non-Water")
    B24_res_idx = groups.index("res-B24")
    B25_res_idx = groups.index("res-B25")
    
    if args.skip >= 1:
        pass
    else:
        os.system(f'mpirun -np 1 gmx_mpi sasa -f {args.xtc} -s {args.tpr} -n {args.ndx} -o {args.sys}_sasa_res_B24.xvg -surface "group {non_water_idx}" -output "group {B24_res_idx}" -dt 250')

    if args.skip >= 2:
        pass
    else:
        os.system(f'mpirun -np 1 gmx_mpi sasa -f {args.xtc} -s {args.tpr} -n {args.ndx} -o {args.sys}_sasa_res_B25.xvg -surface "group {non_water_idx}" -output "group {B25_res_idx}" -dt 250')
