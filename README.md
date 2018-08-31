# PELE_Platform

# PELE_Platform
--------------
Wrap up Platform to launch all PELE features. [AdaptivePELE, MSM, LigandGrowing, Glide Rescoring]

# PELE_Platform Fucntionalities
-------------------------------
1) Calculate absoluite free energies with [MSM_PELE](https://github.com/danielSoler93/MSM_PELE)
2) Run Out-in | Kinase Rescoring | Induce-Fit | Free [Adaptive PELE](https://github.com/AdaptivePELE/AdaptivePELE)
3) Grow your Ligand with [FRAG_PELE] (https://github.com/danielSoler93/LigandGrowing)

# Getting Started
-------------------
0) git clone https://github.com/miniaoshi/PELE_Platform.git


Change all path under PELE_Platform/constants.py

3) Fulfill the next dependencies:
- Adaptive Pele
- pyemma
- msmtools
- Prody 1.8.2
- Pandas

4) Run the platform as:
- python PELE_Platform/main.pdb pdb_with_complex residuename chain --functionality
-  `e.g. python /home/dsoler/PELE_Platform/main.py PR_1A28_xray_-_minimized.pdb STR Z --msm'

# Arguments:
---------------
- **Description:** 
    Monte Carlo Protein Energy Landscape Exploration (PELE) coupled with Markov State Model (MSM) analysis  with the aim to calculate absolute free energies
    - **Requested arguments:** 
        - **complex**: Complex with target & ligand in the binding site.
        - **resname**: Residue name of the ligand in the BS.
        - **chain**: Chain of the ligand in the BS.
    -  `e.g. python /home/dsoler/PelePlop/main.py PR_1A28_xray_-_minimized.pdb STR Z` <br />
    
    - **Optional arguments:**
        - **--charge_ter**       If set charge protein terminals
        - **--forc FORC**        Forcefield to use ["OPLS2005", "AMBER99sb"]
        - **--native NATIVE**    Native file to compare RMSD to
        - **--cpus CPUS**        Number of processors to use in adaptive in out
        - **--core CORE**        Give one atom of the core section for PlopRotTemp
        - **--mtor MTOR**        Gives the maximum number of torsions allowed in each
        group. Will freeze bonds to extend the core if necessary.
        - **--n N**              Maximum Number of Entries in Rotamer File
        - **--clean**            Whether to clean up all PlopRotTemp intermediate files
    - **Output:**
        - The platform will output a file results_summary.txt inside {resname}_PELE/output_pele/results_summary.txt with the           absolute free energy average over all pele trajectories and its standard deviation. If nothing is ouputted means             something is wrong and you can refear to the output.log file or mpi_{IDJOB}.err  mpi_{IDJOB}.out files to traceback           the error.
