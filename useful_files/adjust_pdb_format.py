import argparse 
import os

def initialize():
    parser = argparse.ArgumentParser(
        description="This code adjusts the format of pdb files extracted from MD trajectories to make it GLYCAM-readable.")
    parser.add_argument(
        '-i',
        '--input',
        help='The input PDB file to be reformatted.')
    parser.add_argument(
        '-t',
        '--template',
        help='The GLYCAM-readable PDB template whose coordinates will be replaced by the coordinates of the input PDB file.')

    args_parse = parser.parse_args()

    return args_parse

def check_coordinates(coords_str):
    """
    A sanity check for the parsed string of coordinates.
    """
    coords = coords_str.split()
    for i in coords:
        try:
            float(i)
        except ValueError:
            print("Warning: The coordinates might not be parsed correctly!")


if __name__ == "__main__":
    args = initialize()
    f = open(args.input)
    lines = f.readlines()
    f.close()

    # 1. Get the title for the PDB output (strip() for getting rid of the escape characters)
    title = f'{lines[0].strip()} ({lines[1].split("TITLE     ")[-1].strip()})'

    # 2. Get the coordinates from the input PDB
    coords = []
    n1, n2 = 0, 0  # number of atoms of the two PDB structures
    for l in lines:
        if "ATOM" in l:
            n1 += 1
            coords_str = l[32:54]
            check_coordinates(coords_str)
            coords.append(coords_str)  # A string of coordinates

    # 3. Copy and rename the template and replace the coordinates
    f = open(f"{args.template}")
    lines = f.readlines()
    f.close()

    for i in range(len(lines)):
        if "REMARK" in lines[i]:
            lines[i] = title
        if "ATOM" in lines[i]:
            n2 += 1
            before_coords = lines[i][:32]
            after_coords = lines[i][54:]
            lines[i] = before_coords + coords[n2 - 1] + after_coords

    if n1 != n2:
        print("Warning: Two PDB files might not have the same number of atoms!")
        print(f"The input PDB file {args.input} has {n1} atoms.")
        print(f"The template PDB file {args.template} has {n2} atoms.")

    # save the renewed coordinates
    sys = args.input.split(".pdb")[0]
    with open(f"{sys}_input.pdb", "a") as f:
        for i in range(len(lines)):
            print(lines[i].strip(), file=f)
    f.close()
