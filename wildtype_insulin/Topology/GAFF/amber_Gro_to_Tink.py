import sys
import numpy as np
import itertools as it
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-t', dest='top_file', help='Gromacs topology file')
parser.add_option('-o', dest='output_file', help='String for output .prm file')

(options, args) = parser.parse_args()
top_file = options.top_file
output_file = options.output_file


def extra_tinker_inputs(mass):
    if np.around(mass) == 1.0:
        periodic = '1'
        atom_name = 'H'
    elif np.around(mass) == 12.0:
        periodic = '6'
        atom_name = 'C'
    elif np.around(mass) == 14.0:
        periodic = '7'
        atom_name = 'N'
    elif np.around(mass) == 16.0:
        periodic = '8'
        atom_name = 'O'
    elif np.around(mass) == 32.0:
        periodic = '16'
        atom_name = 'S'
    elif np.around(mass) == 35.0:
        periodic = '17'
        atom_name = 'Cl'
    elif np.around(mass) == 127.0:
        periodic = '53'
        atom_name = 'I'
    else:
        print('ERROR: atom type must be added, exiting code')
        sys.exit()
    return atom_name, periodic


def get_gro_params(params, param):
    found_param = False
    list_params = []
    for i in range(len(params[:, 0])):
        if (params[i, 0] == '[') and (params[i, 1] == param):
            found_param = True
        elif (params[i, 0] == '[') and (params[i, 1] != param):
            found_param = False
        elif found_param:
            if (list(params[i, 0])[0] != ';') and (params[i, 0] != 'nan'):
                list_params.append(i)
    return params[list_params, :]


def max_bonds(bonds_in):
    nbonds = []
    for i in range(np.max(bonds_in[:, :2].astype(int))):
        nbonds.append(len(np.where(i+1 == bonds_in[:, :2].astype(int).flatten())[0]))
    return nbonds


with open(top_file) as f:
    params_in = np.array(list(it.zip_longest(*[lines.split() for lines in f], fillvalue='nan'))).T

# Getting parameters from the topology file
atoms = get_gro_params(params_in, 'atoms')
atom_types = get_gro_params(params_in, 'atomtypes')
bonds = get_gro_params(params_in, 'bonds')
angles = get_gro_params(params_in, 'angles')
dihedrals = get_gro_params(params_in, 'dihedrals')
num_bonds = max_bonds(bonds)

# putting the string together to write out
params_out = "forcefield              AMBER-FF99 \n \nvdwtype                 LENNARD-JONES \n" \
             "radiusrule              ARITHMETIC \nradiustype              R-MIN \n" \
             "radiussize              RADIUS \nepsilonrule             GEOMETRIC \n" \
             "vdw-14-scale            2.0 \nchg-14-scale            1.2 \nelectric                332.0522173" \
             "\ndielectric              1.0\n\n"
s = '  '


# Adding the atoms and atom types
for i in range(len(atoms[:, 0])):
    params_out += "atom" + s + atoms[i, 0] + s + atoms[i, 0] + s
    atom_name, periodic = extra_tinker_inputs(float(atoms[i, 7]))
    params_out += atom_name + s + '"..."' + s + periodic + s + atoms[i, 7] + s + str(num_bonds[i]) + "\n"
params_out += "\n"

# Adding VdW parameters
for i in range(len(atoms[:, 0])):
    params_out += "vdw" + s + atoms[i, 0] + s
    placement = np.where(atoms[i, 1] == atom_types[:, 0])[0][0]
    params_out += str(atom_types[placement, 5].astype(float)*10 / 2 * 2 ** (1/6)) + s
    params_out += str(atom_types[placement, 6].astype(float) / 4.184) + "\n"
params_out += "\n"

# Adding in bonds
for i in range(len(bonds[:, 0])):
    if bonds[i, 2] == '1':
        params_out += "bond" + s + bonds[i, 0] + s + bonds[i, 1] + s
        params_out += str(bonds[i, 4].astype(float) / 4.184 / 10 ** 2 / 2) + s
        params_out += str(bonds[i, 3].astype(float) * 10) + "\n"
    else:
        print("ERROR: bond parameter type " + bonds[i, 2] + " is not supported at the time")
        sys.exit()
params_out += "\n"

# Adding in angles
for i in range(len(angles[:, 0])):
    if angles[i, 3] == '1':
        params_out += "angle" + s + angles[i, 0] + s + angles[i, 1] + s + angles[i, 2] + s
        params_out += str(angles[i, 5].astype(float) / 4.184 / 2) + s
        params_out += str(angles[i, 4].astype(float)) + "\n"
    else:
        print("ERROR: angles parameter type " + angles[i, 3] + " is not supported at the time")
        sys.exit()
params_out += "\n"

# Adding in dihedrals
_, index, inverse = np.unique(dihedrals[:, :5], axis=0, return_index=True, return_inverse=True)
for i in range(len(index)):
    if dihedrals[index[i], 4] == '1':
        params_out += "torsion" + s + dihedrals[index[i], 0] + s + dihedrals[index[i], 1] + s + dihedrals[index[i], 2] \
                      + s + dihedrals[index[i], 3] + s
        for j in np.where(inverse == i)[0]:
            params_out += str(dihedrals[j, 6].astype(float) / 4.184) + s
            params_out += dihedrals[j, 5] + s + dihedrals[j, 7] + s
        params_out += "\n"
    elif dihedrals[index[i], 4] == '4':
        params_out += "imptors" + s + dihedrals[index[i], 0] + s + dihedrals[index[i], 1] + s + dihedrals[index[i], 2] \
                      + s + dihedrals[index[i], 3] + s
        for j in np.where(inverse == i)[0]:
            params_out += str(dihedrals[j, 6].astype(float) / 4.184) + s
            params_out += dihedrals[j, 5] + s + dihedrals[j, 7] + s
        params_out += "\n"
    else:
        print(dihedrals[index[i]])
        print("ERROR: dihedral parameter type " + dihedrals[index[i], 4] + " is not supported at the time")
        sys.exit()
params_out += "\n"

# Adding in charges
for i in range(len(atoms[:,0])):
    params_out += "charge" + s + atoms[i, 0] + s + atoms[i, 6] + "\n"
params_out += "\n"

# Writing Tinker prm file
with open(output_file.split('.')[0] + '.prm', 'w') as file_out:
    file_out.write(params_out)
print("Writing Tinker parameters to: " + output_file.split('.')[0] + '.prm')



















