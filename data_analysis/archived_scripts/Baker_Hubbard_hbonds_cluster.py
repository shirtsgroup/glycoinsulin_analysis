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

class Logging:
    def __init__(self, output):
        self.output = output
    
    def logger(self, *args, **kwargs):
        print(*args, **kwargs)
        with open(self.output, "a") as f:
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


if __name__ == '__main__':
    t1 = time.time()
    args = initialize()
    traj = md.load_pdb(f'glycoform_{args.number}_ACS_cluster_dt250.pdb')
    
    L = Logging(f'ACS_{args.number}_hbonds_cluster_results.txt')
    string = f'Hydrogen bond analysis of glycoform {args.number} cluster medoids'
    L.logger(string)
    L.logger('=' * len(string))
    L.logger(f'Analyzed file: glycoform_{args.number}_ACS_cluster_dt250.pdb')

    ordinal = lambda n: "%d%s" % (n,"tsnrhtdd"[(n//10%10!=1)*(n%10<4)*n%10::4])
    for i in range(5):   # print out the hydrogen bonds for the first five cluster medoids 
        hbond_list = find_glycan_hbond(traj[i], 0)
        if len(hbond_list) == 0:
            L.logger(f'{ordinal(i + 1)} cluster medoid does not have any glycan-involved hydrogen bonds.\n')
        else:
            L.logger(f'{ordinal(i + 1)} cluster medoid has the following glycan-involved hydrogen bonds:\n{",".join(hbond_list)}\n')
