#!/usr/bin/env python3

import argparse
import prody as pd
import xlsxwriter
import os.path
import time

"""
Title of script: dimerOcclusion.py

Author: Dominique A. Ramirez
Version: 2.0
Version Date: 2021 06 25

Description of script:
*This is not a generic script. This is specifically for the glycoinsulin project. This could be extended generally*

This script reads a trajectory file, in PDB format, using ProDy and calculates the number of neighbors within a certain 
distance between the glycan and the dimer interface. It will then spit this out into a Xlsx file for manipulation. It 
will also be set up to read only part of the trajectory if desired.

PLEASE NOTE: the Dimer-Interface residues are hard-coded, as well as the Glycan residues. 

Update logs:
v2.0 - minor updates to clean up the code
"""

""" Set up the ArgParse to get the files that we need """
parser = argparse.ArgumentParser(description="Analyzer of Glycan-Dimer Interface Neighbors from trajectories")

# -- Positional arguments
parser.add_argument("-f", action="store", help="Name of PDB trajectory (in directory) or complete path to file, WITH extension (.pdb)")
parser.add_argument("-o", action="store", help="Name of neighbor Analysis output file with extension (.xlsx)")
parser.add_argument("-d", action="store", help="Distance (in A) for atom neighbor searching")
parser.add_argument("-g", action="store", help="Size of glycan (in residues) i.e. the number of residues that the "
                                               "glycan is assigned to (often 1, 2, or 3)")
parser.add_argument("-t", action="store", help="Time-step size i.e. the time step of each frame (no default provided) (ns)")
parser.add_argument("-fn", action="store", help="(Optional) Every nth frame to analyze, if every frame is not desired")
args = parser.parse_args()

# -- First check that input file is provided, and then check that it exists. Exit if not.
if args.f is None:
    print("No input file given, exiting now.")
    exit(1)
else:
    if os.path.exists(args.f):
        input_file = args.f
    else:
        print("Input file does not exist! (or cannot be found)")
        exit(0)
# -- Second check that the output files are provided. Use default output if none specified. Then check if output already
# -- exists. Exit if so.
if args.o is None:
    output_file = "output.xlsx"
else:
    output_file = args.o
if os.path.exists(output_file):
    print("Output file already exists. Exiting to prevent overwriting data.")
    exit(0)
# -- Third, check if distance is provided and set default if not
if args.d is None:
    distance = 10.0
else:
    try:
        distance = float(args.d)
    except(TypeError):
        print("Distance parameter needs to be a float. Please try again.")
        exit(2)
# -- Fourth, set up the glycan size variable
if args.g is None:
    print("Size of glycan is required. Exiting now.")
    exit(1)
else:
    try:
        glycan_size = int(args.g)
    except(TypeError):
        print("Glycan size parameter must be an integer. Please try again.")
        exit(2)
# -- Fifth, set up the time_step size variable
if args.t is None:
    print("Time-step size is required. Exiting now.")
    exit(1)
else:
    try:
        time_step_size = float(args.t)
    except(TypeError):
        print("Time-step parameter must be a float. Please try again.")
        exit(2)
# -- Lastly, set up the nth frame number
if args.fn is None:
    frame_num = 1
else:
    try:
        frame_num = int(args.fn)
    except(TypeError):
        print("nth Frame to analyze must be an integer. Please try again.")
        exit(2)

start_time = time.time()

""" Perform the distance calculation """
# First parse the PDB file using ProDy, then collect the hierarchical view of the object
traj = pd.parsePDB(input_file)
num_of_frames = traj.numCoordsets()
glyc_hv = traj.getHierView()
#-- Sanity check to make sure that the expected size of the glycan (in residues) matches what is in the file,
# since insulin is 51 residues long.
if (glyc_hv.numResidues() - 51) != glycan_size:
    print("The size of the glycan does not match the number of residues in the trajectory file. Please double "
          "check the size of your glycan and try again.")
    exit(2)

""" Set up the output file before writing data """
workbook = xlsxwriter.Workbook(output_file)
worksheet = workbook.add_worksheet("gly+DI neighbors")

header1 = "Input File Trajectory: " + input_file
header2 = "Distance for neighbor searching: " + str(distance)
worksheet.write_string(0,0,header1)
worksheet.write_string(1,0,header2)

worksheet.write_string(2,0,"Frame No.")
worksheet.write_string(2,1,"Time (ns)")
worksheet.write_string(2,2,"No. Neighbors")
worksheet.write_string(2,3,"Glycan in front of dimer (1 or 0)")
# -- I've already defined content on the first row (0). Since I may not be reading every single frame
# I need to have a separate row counter to keep all my data together in the output file when I go to write things
# below.
row_counter = 3

for i in range(num_of_frames):
    if i % frame_num == 0: #perform modulo to find out if the frame should be analyzed or not
        traj.setACSIndex(i)
        hv = traj.getHierView()
        #-- HARD CODED DIMER INTERFACE
        # This can be changed to any residue combination. Also make sure to specify the chain ID
        # if such a thing exists in the PDB. Format is hv['chainID', residueNumber]
        dimerInterface = hv[' ', 44] + hv[' ', 45] + hv[' ', 46] + hv[' ', 47]

        if glycan_size == 1:
            glycan = hv[' ', 52]
        elif glycan_size == 2:
            glycan = hv[' ', 52] + hv[' ', 53]
        elif glycan_size == 3:
            glycan = hv[' ', 52] + hv[' ', 53] + hv[' ', 54]

        neighbors = pd.findNeighbors(glycan, distance, dimerInterface)
        worksheet.write_string(row_counter,0,str(i))
        worksheet.write_string(row_counter,1,str(round(i*time_step_size,2)))
        worksheet.write_string(row_counter,2,str(len(neighbors)))
        if len(neighbors) > 0:
            worksheet.write_string(row_counter,3,str("1"))
        else:
            worksheet.write_string(row_counter,3,str("0"))

        row_counter += 1
    else:
        continue

workbook.close()

elapsed_time = time.time() - start_time
<<<<<<< HEAD
print(time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
=======
print(time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
>>>>>>> d298c0e3477e4cadd6ddf9852d9e8a6d49cb7ed8
