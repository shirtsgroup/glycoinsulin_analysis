In the folder `Glycam_outputs`, the final outputs of GLYCAM, `structure.parm7` and `structure.rst7` were converted by ACPYPE to `structure_GMX.top` and `structure_GMX.gro` by the following command:
```
acpype -p structure.parm7 -x structure.rst7
```
In the folder `Topology`, there are `glycoform_14.gro` and `glycoform_14.top`. They are just copies of `structure_GMX.gro` and `structure_GMX.top` and are the inputs of the procedure starting from the folder `Box`.



