# Protocol for preparing input files for MD simulations
Below we summarize the protocol we used to prepare the input files for the MD simulations of insulin wildtypes and glycoforms. Note that typically the first two steps are finished locally, while the remaining steps are carried out on an HPC after file transferring.  
- **Step 1**: Parameterize the structure of the insulin wildtype which using H++
- **Step 2**: Convert the file formats
  - *For insulin wildtype simulations*:
    Simply convert the AMBER coordinates and topology files output by H++ into GROMACS formats using `acpype.py`:
    ```
    python acpype.py -x {sys}_Hpp.crd -p {sys}_Hpp.top
    ```
  - *For insulin glycoform simulations*:
  Use GLYCAM glycoprotein builder to build glycans attached to the insulin wildtype. Download the output files and convert the AMBER coordinates and topology files into GROMACS formats using `acpype`:
    ```
    python acpype.py -x structure.rst7 -p structure.parm7
    ``` 
- **Step 3**: Execute the following command
  ```
  bash prep_dir.sh
  ```
- **Step 4**: Modify the topology file. 
  - Add the following lines right before the `[ system ]` directive:
    ```
    ; Include water topology
    #include "amber99sb-ildn.ff/tip3p.itp"

    ; Include topology for ions
    #include "amber99sb-ildn.ff/ions.itp"
    ```
  - Append the following atom types:
    ```
    OW       OW          16.00    0.0000    A     3.15061e-01  6.36386e-01
    HW       HW          1.008    0.0000    A     0.00000e+00  0.00000e+00
    Cl       Cl          35.45    0.00000   A    4.40104e-01    4.18400e-01
    Na       Na          22.99    0.00000   A    3.32840e-01    1.15897e-02
    IB       IB          0.00000  0.00000   A    8.90899e-01    4.18400e-01
    C0       C0          0.00000  0.00000   A    3.05240e-01    1.92376e+00
    MG       MG          0.00000  0.00000   A    1.41225e-01    3.74342e+00
    K        K           0.00000  0.00000   A    4.73602e-01    1.37235e-03
    Rb       Rb          0.00000  0.00000   A    5.26699e-01    7.11280e-04
    Cs       Cs          0.00000  0.00000   A    6.04920e-01    3.37230e-04
    Li       Li          0.00000  0.00000   A    2.02590e-01    7.65672e-02
    Zn       Zn          0.00000  0.00000   A    1.95998e-01    5.23000e-02
    ```
- **Step 5**: On Bridges-2, launch an interactive node using: `interact -N 1 --ntasks-per-node=64` and then execute `bash prep_md.sh`.
- **Step 6**: Submit the job as needed.

For more information, please visit the relevant [HackMD note](https://hackmd.io/@WeiTseHsu/glycoinsulin_preparation).






