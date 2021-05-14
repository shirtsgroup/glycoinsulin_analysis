# Protocol for preparing input files for MD simulations
Below we summarize the protocol we used to prepare the input files for the MD simulations of insulin wildtypes and glycoforms. Note that typically the first two steps are finished locally, while the remaining steps are carried out on an HPC after file transferring.  
- **Step 1**: Parameterize the structure of the insulin wildtype of interest or the wildtype which the glycoform of interest is based on using H++.
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
- **Step 3**: Transfer the following files to the working directory on Bridges-2:
  - `mdp_files/`
  - `gmx_analysis.sh`
  - `prep.sh`
- **Step 3**: On Bridges-2, launch an interactive node using: `interact -N 1 --ntasks-per-node=64` and then execute `bash prep.sh`, answer all the prompts to launch an MD simulation, and submit the job. 

For more information, please visit the relevant [HackMD note](https://hackmd.io/@WeiTseHsu/glycoinsulin_preparation).






