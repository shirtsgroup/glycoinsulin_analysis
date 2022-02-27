import numpy as np
import MDAnalysis as mda 

def get_charges(top):
    """
    This function gets the following two lists from a topology:
    (1) A list of total charges up to the n-th atom, where n ranges from 1 to the total numbere of atoms.
    (2) A list of the charge of each atom
    """
    infile = open(top, 'r')
    lines = infile.readlines()
    infile.close()

    qtot = None
    qtot_list = []
    for l in lines:
        if 'qtot' in l and 'nr' not in l:
            qtot = float(l.split(';')[1].split('qtot')[-1])
            qtot_list.append(qtot)
        if qtot is not None and 'qtot' not in l:
            break
    qtot = np.array(qtot_list)
    q = np.insert(np.diff(qtot), 0, qtot[0])

    return qtot, q

def compare_structures(gro_1, gro_2, top_1, top_2, name_1, name_2):
    u_ref = mda.Universe(gro_1)
    u = mda.Universe(gro_2)
    res_ref = u_ref.residues
    res = u.residues

    s_ref = u_ref.select_atoms("all")
    s = u.select_atoms("all")
    qtot_ref, q_ref = get_charges(top_1)
    qtot, q = get_charges(top_2)

    ordinal = lambda n: "%d%s" % (n,"tsnrhtdd"[(n//10%10!=1)*(n%10<4)*n%10::4])
    print(f'-- [ Comparison between {name_1} and {name_2}] --')

    # Step 1: Check the number of atoms
    if len(s_ref) != len(s):
        print(f'{name_1} and {name_2} have different numbers of atoms ({len(s_ref)} v.s. {len(s)})!')
        if len(res_ref) != len(res):
            print(f'    {name_1} and {name_2} have different numbers of residues!')
        else:
            for i in range(len(res_ref)):
                len_1 = len(u_ref.select_atoms(f"resid {i + 1}"))
                len_2 = len(u.select_atoms(f"resid {i + 1}"))
                if len_1 != len_2:
                    # Note that resid starts from 1
                    print(f'    {name_1} and {name_2} have different numbers of atoms at the {ordinal(i + 1)} residue:')
                    print(f'    {name_1}: {len_1} atoms at {res_ref[i]} / {name_2}: {len_2} atoms at {res[i]}')
        print(f'--> \U0000274E Conclusion: {name_1} and {name_2} have different structures!\n')
    else:
        print(f'{name_1} and {name_2} have the same numbers of atoms ({len(s_ref)})!')
        
        # Step 2: Check the total charges
        if qtot_ref[-1] == qtot[-1]:
            print(f'{name_1} and {name_2} have the same total charges ({qtot_ref[-1]})!')
        else:
            print(f'{name_1} and {name_2} have different total charges!')

        # Step 3: Check the atom type for each atom
        score = 0
        resid_diff = []
        res_diff = None   # check if any reisdues are different between structures
        for i in range(len(res_ref)):
            if len(res_ref[i].atoms) != len(res[i].atoms):  # check if # of atoms are the same for residue i
                res_diff = True  

            n = np.min([len(res_ref[i].atoms), len(res[i].atoms)])  # might have different # of atoms
            for j in range(n):
                score += (res_ref[i].atoms[j].name == res[i].atoms[j].name)
                if (res_ref[i].atoms[j].name == res[i].atoms[j].name) is False:
                    if i not in resid_diff:
                        resid_diff.append(i)
                else:  
                    # Atoms of the same atom types in differently protonated residues would have different charges
                    # Here we compare the atoms with the same atom types and residue types 
                    idx_1, idx_2 = res_ref[i].atoms[j].ix, res[i].atoms[j].ix
                    if res_ref[i].resname == res[i].resname: 
                        if np.isclose(q_ref[idx_1], q[idx_2], rtol=1e-05, atol=1e-08) is False:
                            print('    Atoms with the same atom and residue types had different charges!')

        if score == len(s_ref) and res_diff is None:
            print(f'{name_1} and {name_2} have exactly the same atom type for each atom.')
            print(f'--> \U00002705 Conclusion: {name_1} and {name_2} have the same structures.\n')
        else:
            print(f'{name_1} and {name_2} do not have the same atom types.')
            print('    Below is a list of residues that had different atom types:')
            for r in resid_diff:
                print(f'    {name_1}: {res_ref[r]} / {name_2}: {res[r]}')
            print(f'--> \U0000274E Conclusion: {name_1} and {name_2} have different structures!\n')

if __name__ == "__main__":
    WTs = ['4EYD', '4EY9', '4EY1', '3I3Z', '2MVC']
    # Section 1: Compare the wildtype structures
    string = 'Section 1: Comparison of the wildtype structures'
    print(string)
    print('=' * len(string))
    gro_1 = '../../wildtype_insulin/4EYD/Hpp_results/4eyd.gro'
    top_1 = '../../wildtype_insulin/4EYD/Hpp_results/4eyd.top'
    for i in range(1, 5):
        gro_2 = f'../../wildtype_insulin/{WTs[i]}/Hpp_results/{WTs[i].lower()}.gro'
        top_2 = f'../../wildtype_insulin/{WTs[i]}/Hpp_results/{WTs[i].lower()}.top'
        compare_structures(gro_1, gro_2, top_1, top_2, 'WT/4EYD', f'WT/{WTs[i]}')

    # Section 2: Compare the insulin glycoforms
    string = '\nSection 2: Comparison between insulin glycoforms'
    print(string)
    print('=' * len(string))
    for i in range(2, 14):
        gro_1 = f'4EYD_glycoforms/glycoform_{i}_ACS/Glycam_outputs/structure_GMX.gro'
        top_1 = f'4EYD_glycoforms/glycoform_{i}_ACS/Glycam_outputs/structure_GMX.top'
        for j in range(1, 5):
            gro_2 = f'{WTs[j]}_glycoforms/glycoform_{i}_ACS/Glycam_outputs/structure_GMX.gro'
            top_2 = f'{WTs[j]}_glycoforms/glycoform_{i}_ACS/Glycam_outputs/structure_GMX.top'
            compare_structures(gro_1, gro_2, top_1, top_2, f'GF{i}/4EYD', f'GF{i}/{WTs[j]}')
