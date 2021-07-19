# Glycoforms of insulin
Below is the directory structure of the folder `glyco_insulin`, with `i` being the index number of the glycoform:
```
.
├── 2MVC_glycoforms
│   ├── 0.15_80_10_pH7.2_2mvc_clean_noH.result.pdb
│   ├── 2mvc_input.pdb
│   └── glycoform_${i}_ACS
├── 3I3Z_glycoforms
│   ├── 0.15_80_10_pH6.9_3i3z_clean_mutated.result.pdb
│   ├── 3i3z_input.pdb
│   └── glycoform_${i}_ACS
├── 4EY1_glycoforms
│   ├── 0.15_80_10_pH7.9_4ey1_clean.result.pdb
│   ├── 4ey1_input.pdb
│   └── glycoform_${i}_ACS
├── 4EY9_glycoforms
│   ├── 0.15_80_10_pH8.0_4ey9_clean.result.pdb
│   ├── 4ey9_input.pdb
│   ├── 4eyd_max_sasa_SerA12_input.pdb
│   └── glycoform_${i}_ACS
├── 4EYD_glycoforms
│   ├── 0.15_80_10_pH8.0_4eyd_clean.result.pdb
│   ├── 4eyd_input.pdb
│   ├── 4eyd_max_sasa_SerA12_input.pdb
│   ├── archived
│   └── glycoform_${i}_ACS
├── README.md
├── compare_results.txt
└── compare_structures.py
```
Each folder `{PDB_ID}_glycoforms` stores two PDB files and a series of folder for different glycoforms. The PDB file with the longer file name (in the form of `{salinty}_{ext_dielectric_const}_{int_dielectric_const}_pH{pH_value}_{PDB_ID}_clean.result.pdb`) is the PDB file output from the H++ server under the desired conditions (mainly the protonation states), while the other one (file name in the form of `{PDB_ID}_input.pdb`) is the PDB file used as the input to GLYCAM glycoprotein builder. The only difference between the two PDB file is that the second one has a chain identifier for each residue, which is recommended when attaching sugar molecules.

As shown in the directory structure above, various glycoforms are based on 5 different X-ray or NMR-resolved structures, with PDB code as 2MVC, 3I3Z, 4EY1, 4EY9, and 4EYD. For each of them, we prepared 12 glycoforms as shown in the figure below, which were the systems of interest in the [2018 ACS paper by Tan et al](https://pubs.acs.org/doi/abs/10.1021/acschembio.7b00794).

<center><img src=https://i.imgur.com/kwdYFCt.png width=400>
</center>

In each subfolder `glycoform_{n}_ACS`, there are two PDB files. One is the raw output from the H++ server, while the other one was the PDB file we directly use as the input to the GLYCAM glycoprotein builder. To prepare the latter, we manually edited the former to assign chain A and chain B and strictly followed the conventional [PDB format](https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html).

For each glycoform, below we summarize the details of parameters adopted in the parameterization process done by the GLYCAM glycoprotein builder. Note that for structure 4EYD, and 4EY9, the glycosylation site Ser A12 was predicted to be buried in the protein and GLYCAM was not able to attach the glycan. (GLYCAM concluded that it is unlikely to have glycosylation there.) Therefore, for glycoform 3, instead of directly using the PDB file output by H++, we drew the configuration whose SASA was the largest from the MD simulation of the wildtype to use as the input to GLYCAM. 
### 1. Glycoforms 2 to 6 in the 2018 ACS paper
- Attached glycan: all GalNAc ($\alpha$ form, D-form, with OH as the aglycon)
- Glycosylation sites: SerA9, ~~SerA12~~, SerB9, ThrB27, ThrB30
- No glycosidic linkage
- GLYCAM notation: all `DGalpNAca1-OH`

### 2. Glycoforms 7, 9, and 11 in the 2018 ACS paper
- Attached glycan: all Man ($\alpha$-form, D-form, with OH as the aglycon)
- Glycosylation sites: SerA9, ThrB27, ThrB30
- No glycosidic linkage
- GLYCAM notation: all `DManpa1-OH`

### 3. Glycoforms 8, 10, and 12 in the 2018 ACS paper 
- Attached glycan: all di-mannose ($\alpha$-form, D-form, with OH as the aglycon)
- Glycosylation sites: SerA9, ThrB27, ThrB30
- Linkage: all $\alpha$-1, 2-glycosidic linkage
- Default glycosidic angles were used.
- GLYCAM notation: all `DManpa1-2DManpa1-OH`

### 4. Glycoforms 13 in the 2018 ACS paper
- Manα2Manα2Manα-ThrB27
- Attached glycan: a tri-mannose mannose ($\alpha$-form, D-form, with OH as the aglycon)
- Glycosylation site: ThrB27
- Linkage: all $\alpha$-1, 2-glycosidic linkage
- Default glycosidic angles were used.
- GLYCAM notation: `DManpa1-2DManpa1-2DManpa1-OH`

**Note**: For the X-ray crystal strucutres of the wildtype 4EYD and 4EY9, the solvent accessible surface area of the residue SerA12 was too small for GLYCAM glycoprotein builder to attach any sugar. Therefore, from the MD trajectories of these two wildtype strucutres, we extracted the configuration with the largest SASA value at SerA12. The extracted configuration went through the some postprocessing to serve as the input for the GLYCAM glycoprotein builder to attach GalNAc at SerA12. The detailed workflow is as follows:
- For each wildtype structure of 4EYD and 4EY9, make the index group for the residue SerA12.
  - Sample command: `gmx make_ndx -f ../MD/*gro -o *ndx -n *ndx`
- Calculate the SASA of SerA12
  - Sample command: `gmx sasa -f ../MD/traj_comp.xtc -s ../MD/4ey9_md.tpr -n 4ey9.ndx -o 4ey9_sasa_SerA12.xvg -surface 'group 14' -output 'group 21' -dt 250`
- Find the time frame corresponding to the maximum of the SASA.
- Use `trjconv` to recenter the protein to avoid crossing the periodic boundaries. Note that this is necessary since GLYCAM glycoprotein builder does not seem to take PBC into account. If important structures (e.g. a disulfide bridge) are crossing the periodic boundaries, glycoprotein builder would not be able to catch it.
  - Sample command 1: `gmx trjconv -s ../MD/*tpr -f ../MD/traj_comp.xtc -o 4eyd_nojump_protein_dt2.xtc -center -pbc nojump`
  - Sample command 2: `gmx trjconv -s ../MD/*tpr -f *nojump*.xtc -o 4eyd_center_protein_dt2.xtc -center -pbc mol -ur compact`
- Use `trjconv` to extract the configuration PDB file (output: only protein)
  - Sample command: `gmx trjconv -s ../MD/*tpr -f 4eyd_center_protein_dt2.xtc -o 4eyd_max_sasa_SerA12.pdb -dump {time}` where `{time}` is the timeframe of the configuration of interest corresponds.
- Modify the output PDB file to make it readable by GLYCAM glycoprotein builder
  - The easiest way is to copy the input file for GLYCAM glycoprotein builder when attaching other sugars to other sites, i.e. `4eyd_input.pdb` and `4ey9_input.pdb` and then replace the coordinates with the coordinates of the extracted configuration. 
  - Once can modify `REMARK` to make the source of the PDB file more clear. However, statements like `TITLE` or `CRYSTL1`, `MODEL` are not compatible in GLYCAM.
- With inputs such as `4eyd_max_sasa_SerA12_input.pdb` or `4ey9_max_sasa_SerA12_input.pdb`, attach GalNAc at SerA12.
- Similarly, use `prep.sh` to automate the preparation process and launch the MD simulation. 