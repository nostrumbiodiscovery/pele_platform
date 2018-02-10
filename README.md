# MSM_Pele
--------------
Monte Carlo Protein Energy Landscape Exploration (PELE) coupled with Markov State Model (MSM) analysis  with the aim to calculate absolute free energies.

# MSM_Pele's PipeLine
-------------------------------
1) [Protein Preparation for Pele](https://github.com/Jelisa/mut-prep4pele)
2) [PlopRotTemp_SCHR2017](https://github.com/miniaoshi/PlopRotTemp_S_2017)
3) [Adaptive PELE](https://github.com/AdaptivePELE/AdaptivePELE)
4) [PELE(comercial software)](https://pele.bsc.es/pele.wt)
5) [MSM](https://github.com/miniaoshi/Pele_scripts)

# Getting Started
-------------------
0) git clone https://github.com/miniaoshi/MSM_PELE.git
1) Change all paths under **MSM_Pele/constants.py** to your local environment.

2) Change LD_LYBRARY_PATH as next in the conf.sl to make schrodinger python work:
-  `e.g. export LD_LIBRARY_PATH=**/sNow/easybuild/centos/7.4.1708/Skylake/software/schrodinger2017-4/mmshare-v4.0/lib/Linux-x86_64/**:$LD_LIBRARY_PATH`
3) Fulfill the next dependencies:
- Adaptive Pele
- pyemma
- msmtools
- Prody 1.8.2
- Pandas
4) Run the platform as:
- python MSM_PELE/main.pdb pdb_with_complex residuename chain
-  `e.g. python /home/dsoler/PelePlop/main.py PR_1A28_xray_-_minimized.pdb STR Z`
- slurm: `sbatch conf.sl` (change the complex path)

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
