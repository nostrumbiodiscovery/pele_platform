#!/bin/bash

#Parse arguments from input
source $(dirname $0)/Helpers/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('pdb_file', type=str, help="PDB of COMPLEX to run PELE on.")
parser.add_argument('ligand_chain', type=str, help="Cain of the ligand to be parametrizeo on the pdb_file")
EOF
if [ "$SCHRODINGER" == "" ]; then
	echo "SCHRODINGER IS NOT EXPORTED"
else

	#Get ligand name from mae input file
	pdb=$(basename "$PDB_FILE")
	pdbname="${pdb%.*}"

    #Global varibales
    ligand_mae="${pdbname}.mae"

    #extract ligand
	mapfile -t ligand_pdb < <(python Helpers/extract_ligand.py --pdb $pdb --general_name $pdbname --ligand_chain $LIGAND_CHAIN --executing_folder $PWD)

	#PlopRotTemp over ligand
	$SCHRODINGER/utilities/structconvert -ipdb ${ligand_pdb[0]} -omae $ligand_mae
	$SCHRODINGER/utilities/python PlpRotTemp/main.py $ligand_mae
	rm ${ligand_pdb[0]}
	mapfile -t output_files < input.txt

	#Prepare PELE env
    Pele_directory="${pdbname}_Pele"
    echo "$Pele_directory"
    if [ ! -d "$Pele_directory" ]; then
		mkdir $Pele_directory
	fi
	
	if [ ! -d "${Pele_directory}/DataLocal/Templates/OPLS2005/HeteroAtoms/" ]; then
		mkdir -p "${Pele_directory}/DataLocal/Templates/OPLS2005/HeteroAtoms/"
	fi
	mv "${output_files[0]}" "${Pele_directory}/DataLocal/Templates/OPLS2005/HeteroAtoms/"
	
	if [ ! -d "${Pele_directory}/DataLocal/LigandRotamerLibs" ]; then
		mkdir -p "${Pele_directory}/DataLocal/LigandRotamerLibs"
	fi
	mv "${output_files[1]}" "${Pele_directory}/DataLocal/LigandRotamerLibs"

	if [ ! -d "${Pele_directory}/results" ]; then
		mkdir "${Pele_directory}/results"
	fi

	cp $pdb "${Pele_directory}/complex.pdb"
	cp $pdb "${Pele_directory}/native.pdb"

	cp PeleTemplates/*control_file* $Pele_directory

	rm ligand_mae
	rm input.txt

	cd "$Pele_directory"
	ln -s /home/dani/repos/PELErev12535/Data/ Data
	ln -s /home/dani/repos/PELErev12535/Documents/ Documents

	#RunPele
	/usr/lib64/openmpi/bin/mpirun -np 3 /opt/PELE/v12354/bin/PELE-1.5_mpi control_file --license-directory ~/repos/PELErev12535/licenses/

fi








