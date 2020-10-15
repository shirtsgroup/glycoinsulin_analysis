## How to parameterize an arbitrary monomer with GAFF using antechamber ##
# Updated Methodology #
1. Create structure in ChemDraw (or any other tool i.e. Avogadro). Save as *.mol (choose MDL molfile, not MDL molfile V3000)
2. Convert to .pdb using <a href="https://cactus.nci.nih.gov/translate/">online smiles converter</a> or openbabel.
3. Make sure the .pdb has the residue name in the correct location (starts at character 18 in ATOM section). 
    1. If it is blank, you will get a difficult-to-trace memory allocation error like "_memory allocation error for *bond_array_"
4. Run GAFF_params.py: python GAFF_params.py -p *.pdb -o output -m *.mdp -n <netcharge>
    1. Provide the flag --tinker to then convert the Gromacs topolgy file to a Tinker topology file.

_NOTE: When converting from Gromacs to Tinker verify that the energies match up between programs. One common error is that Tinker does not get the correct number of improper torsions and the atom order in the .prm file needs to be reversed._

# Old Methodology #
Follow steps 1 -- 3 in the Updated Methodology section.

4. Parameterize with antechamber following the <a href="http://ambermd.org/tutorials/basic/tutorial4b/">tutorial</a>.
5. Convert the output *.prmtop and *.inpcrd files to gromacs *.gro and *.top file using acpype.py
6. Energy minimize the structure
7. Reassign partial charges using OpenEye software
    1. Use <a href="https://docs.eyesopen.com/quacpac/molchargeusage.html">molcharge</a> which is shipped with QUACPAC
    2. Use the <a href="https://docs.eyesopen.com/quacpac/molchargetheory.html#am1bcc-charges">am1bccsym</a> method to assign charges symmetrically after a conformational search is performed
    3. Use insertmol2charges.py to replace the charges output by antechamber with those output by molcharge
8. Energy minimize again to get final structure with charges and parameters

_NOTE: Steps 4 - 8 can be done automatically using param.sh_
