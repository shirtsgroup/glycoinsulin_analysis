import time
import numpy as np
import argparse
import mdtraj as md


def initialize():
    parser = argparse.ArgumentParser(
        description='This code uses Baker-Hubbard criterion to identify H.B. for a speicifed glycoform given its trajectory in as a pdb file.')
    parser.add_argument('-n',
                        '--number',
                        type=int,
                        help='The index of the glycoform.')
    args_parse = parser.parse_args()

    return args_parse


def logger(output, *args, **kwargs):
    print(*args, **kwargs)
    with open(output, "a") as f:
        print(file=f, *args, **kwargs)


def find_glycan_hbond(t, freq):
    """
    Parameters
    ----------
    t    (mdtraj Trajectory object) : an mdtraj object generated after loading 
         the trajectory by using md.load_pdb().
    freq (float) : Return only hydrogen bonds that occur in greater this 
         fraction of the frames in the trajectory.

    Return
    ------
    glycan_hbonds (list): A list of string showing donor--acceptor 
                  pairs of all glycan-involved H.B.'s.
    """
    hbonds_idx = md.baker_hubbard(t, freq=freq, periodic=False)
    def convert_idx(
        hbond): return '%s -- %s' % (t.topology.atom(hbond[0]), t.topology.atom(hbond[2]))
    hbonds = []
    for i in hbonds_idx:
        hbonds.append(convert_idx(i))

    glycan_hbonds = []
    for i in hbonds:
        if int(i.split(' -- ')[1].split('-')[0][3:]) > 51:  # glycan as the acceptor
            glycan_hbonds.append(i)

    return glycan_hbonds


def glycan_hbonds_info(t):
    """
    This function estimates the fraction of the occurance of glcan-involved H.B.'s 
    given a mdtraj Trajedtory object.

    Parameter
    ---------
    t    (mdtraj Trajectory object) : an mdtraj object generated after loading 
         the trajectory by using md.load_pdb().

    Return
    ------
    hbond_info (dict): a dictionary with keys as D-A pairs and values as 
               corrsponding fractions of occurence.
    """
    hbond_info = {}
    freq = np.arange(0.10, 1.01, 0.01)
    for i in range(len(freq)):
        hbond_low = set(find_glycan_hbond(t, freq[i]))
        hbond_high = set(find_glycan_hbond(t, freq[i + 1]))
        if len(hbond_high) < len(hbond_low):
            complement_set = hbond_low - hbond_high
            for j in complement_set:
                hbond_info[j] = f'{freq[i]: .2f}'
        if len(hbond_high) == 0:
            break

    result = {k: float(v) for k, v in sorted(
        hbond_info.items(), key=lambda item: item[1], reverse=True)}

    return result


if __name__ == '__main__':
    t1 = time.time()
    args = initialize()
    traj = f'glycoform_{args.number}_ACS_dt250.pdb'
    t = md.load_pdb(traj)
    result = glycan_hbonds_info(t)
    string = f'Hydrogen bond analysis of glycoform {args.number}'
    logger(f'ACS_{args.number}_hbonds_results.txt', string)
    logger(f'ACS_{args.number}_hbonds_results.txt', '=' * len(string))
    logger(f'ACS_{args.number}_hbonds_results.txt', f'Analyzed file: {traj}')

    if len(result) == 0:
        logger(f'ACS_{args.number}_hbonds_results.txt',
               'No glycan-involved hydrogen bonds existed more than 10% of the trajectory.')
    else:
        logger(f'ACS_{args.number}_hbonds_results.txt',
               'Glycan-involved hydrogen bonds:')
        for i in result:
            logger(f'ACS_{args.number}_hbonds_results.txt',
                   f'{i} ({int(result[i] * 100)}%)')
    logger(f'ACS_{args.number}_hbonds_results.txt',
           f'Elapsed time: {time.time() - t1} seconds.')
