import numpy as np
import MDAnalysis as mda 

def get_total_charges(top):
    """
    This function gets the total charges up to the n-th atom from 
    a topology, where n ranges from 1 to the total numbere of atoms.
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
    qtot_list = np.array(qtot_list)

    return qtot_list

def compare_structures(gro_1, gro_2, top_1, top_2, name_1, name_2):
    u_ref = mda.Universe(gro_1)
    u = mda.Universe(gro_2)
    res_ref = u_ref.residues
    res = u.residues


    s_ref = u_ref.select_atoms("all")
    s = u.select_atoms("all")
    qtot_ref = get_total_charges(top_1)
    qtot = get_total_charges(top_2)

    print(f'-- [ Comparison between {name_1} and {name_2}] --')
    if len(s_ref) != len(s):
        print(f'{name_1} and {name_2} have different numbers of atoms!')
        print(f'--> Conclusion: {name_1} and {name_2} have different strucutres!\n')
    else:
        print(f'{name_1} and {name_2} have the same numbers of atoms!')
        if qtot_ref[-1] == qtot[-1]:
            print(f'{name_1} and {name_2} have the same total charges!')
        else:
            print(f'{name_1} and {name_2} have different total charges!')

        score = 0
        resid_diff = []
        res_diff = None   # check if any reisdues are different between structures
        for i in range(len(res_ref)):
            if len(res_ref[i].atoms) != len(res[i].atoms):
                res_diff = True     
            n = np.min([len(res_ref[i].atoms), len(res[i].atoms)])  # might have different # of atoms
            for j in range(n):
                score += (res_ref[i].atoms[j].name == res[i].atoms[j].name)
                if (res_ref[i].atoms[j].name == res[i].atoms[j].name) is False:
                    if i not in resid_diff:
                        resid_diff.append(i)
        if score == len(s_ref) and res_diff is None:
            print(f'{name_1} and {name_2} have exactly the same atom type for each atom.')
            print(f'--> Conclusion: {name_1} and {name_2} have the same structures.\n')
        else:
            print(f'{name_1} and {name_2} do not have the same atom types.')
            print('Below is a list of residues that had different atom types:')
            for r in resid_diff:
                print(f'{name_1}: {res_ref[r]} / {name_2}: {res[r]}')
            print(f'--> Conclusion: {name_1} and {name_2} have different structures!\n')

if __name__ == "__main__":
    WTs = ['4EYD', '4EY9', '4EY1', '3I3Z', '2MVC']
    # Section 1: Compare the wildtype strucutres
    string = 'Section 1: Comparison of the wildtype structures'
    print(string)
    print('=' * len(string))
    gro_1 = '../wildtype_insulin/4EYD/Hpp_results/4eyd.gro'
    top_1 = '../wildtype_insulin/4EYD/Hpp_results/4eyd.top'
    for i in range(1, 5):
        gro_2 = f'../wildtype_insulin/{WTs[i]}/Hpp_results/{WTs[i].lower()}.gro'
        top_2 = f'../wildtype_insulin/{WTs[i]}/Hpp_results/{WTs[i].lower()}.top'
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
        