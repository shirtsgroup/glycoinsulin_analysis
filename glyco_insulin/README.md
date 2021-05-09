# Glycoforms of insulin
Here we store the input files (output by GLYCAM with some postprocessing) of the simulations of various glycoforms based on 5 different X-ray or NMR-resolved structures, with PDB code as 2MVC, 3I3Z, 4EY1, 4EY9, and 4EYD. For each of them, we prepared 12 glycoforms as shown in the figure below, which were the systems of interest in the [2018 ACS paper by Tan et al](https://pubs.acs.org/doi/abs/10.1021/acschembio.7b00794).

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