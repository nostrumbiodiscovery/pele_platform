#!/bin/bash

#######################################################
# PelePlop: Protein Energy Landscape Exploration Platform
# 
# https://github.com/miniaoshi/PelePlop
#
# Author: Daniel Soler, daniel.soler@nostrumbiodiscovery.com
#
# Copyright 2017 BSC, NBD
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
########################################################

PlopRotTemp=$(dirname $0)
FORCEFIELD="OPLS2005"

####################### Check env variables are set ######################

if [[ $(which mpirun 2>&1 > /dev/null) != "" ]]; then
	echo "set mpirun $PATH with: $: set export PATH=/path/to/binary/:$PATH"
	exit 1
elif [[ $(which PELE-1.5_mpi 2>&1 > /dev/null) != "" ]]; then
	echo "set PELE-1.5 binary folder to $PATH with: $: set export PATH=/path/to/binary/:$PATH"
	exit 1
fi

###################### Parse arguments from input ######################

source "${PlopRotTemp}"/Helpers/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('pdb_file', type=str, help="PDB of COMPLEX to run PELE on.")
parser.add_argument('ligand_residue', type=str, help="Residue of the ligand to be parametrizeo on the pdb_file")
parser.add_argument('ligand_chain', type=str, help="Chain of the ligand to be parametrizeo on the pdb_file")
parser.add_argument("--mtor", type=int, help="Gives the maximum number of torsions allowed in each \
                     group.  Will freeze bonds to extend the core if \
                     necessary.")
parser.add_argument("--core", type=int, help="Give one atom of the core section")
parser.add_argument("--n", type=int, help="Maximum Number of Entries in Rotamer File")
parser.add_argument("--mae_charges", help="Use charges specified in the ligand.mae file", action='store_true')
parser.add_argument("--clean", help="To clean up all the intermediate files", action='store_true')
parser.add_argument("--cpus", type=int, help="Number of cores the progam will try to use")
parser.add_argument("--confile", type=str, help="Your own PELE conf file")
parser.add_argument("--native", type=str, help="Native structure to calculate the RMSD")
parser.add_argument("--forcefield", type=str, help="Forcefield to be used. (default=OPLS2005, AMBER99, AMBERBSC)")
EOF

PlopRotTemp_opt=""
Pele_opt=""
if [ "$MTOR" != "" ]; then
	PlopRotTemp_opt="$PlopRotTemp_opt --mtor ${MTOR} "
fi

if [ "$CORE" != "" ]; then
	PlopRotTemp_opt="$PlopRotTemp_opt --core ${CORE}"
fi

if [ "$N" != "" ]; then
	PlopRotTemp_opt="$PlopRotTemp_opt --n ${N}"
fi

if [ "$MAE_CHARGES" != "" ]; then
	PlopRotTemp_opt="$PlopRotTemp_opt --mae_charges"
fi

if [ "$CLEAN" != "" ]; then
	PlopRotTemp_opt="$PlopRotTemp_opt --clean"
fi

if [ "$PlopRotTemp_opt" != "" ]; then
	echo "Using PlopRotTemp options ${PlopRotTemp_opt}"
fi

if [ "$CPUS" != "" ]; then
	echo "Using ${CPUS} cpus"
else
	CPUS=3
fi

if [ "$CONFILE" != "" ]; then
	echo "Using control file $CONFILE"
fi

if [ "$NATIVE" != "" ]; then
	echo "Using native structure $NATIVE"
else
	NATIVE="$PDB_FILE"
	echo "Using input file as the native structure to compute RMSD"
fi
####################################################################################


############################### RUN PELEPLOP ########################################

if [ "$SCHRODINGER" == "" ]; then
	echo "SCHRODINGER IS NOT EXPORTED"
else

	#Get ligand name from mae input file
	pdb=$(basename "$PDB_FILE")
	pdbname="${pdb%.*}"

    #Global varibales
    ligand_mae="${pdbname}.mae"

    #extract ligand
	mapfile -t ligand_pdb < <(python "${PlopRotTemp}/Helpers/extract_ligand.py" --pdb $PDB_FILE --general_name $pdbname --ligand_chain "${LIGAND_RESIDUE} ${LIGAND_CHAIN}" --executing_folder $PWD)

	#PlopRotTemp over ligan
	$SCHRODINGER/utilities/structconvert -ipdb ${ligand_pdb[0]} -omae $ligand_mae
	echo "$SCHRODINGER/utilities/python "${PlopRotTemp}/PlpRotTemp/main.py" $PlopRotTemp_opt $ligand_mae"
	$SCHRODINGER/utilities/python "${PlopRotTemp}/PlpRotTemp/main.py" $PlopRotTemp_opt $ligand_mae
	rm ${ligand_pdb[0]}
	mapfile -t output_files < input.txt

	#Prepare PELE env
    Pele_directory="${pdbname}_Pele"
  
    if [ ! -d "$Pele_directory" ]; then
		mkdir $Pele_directory
	fi
	
	if [ ! -d "${Pele_directory}/DataLocal/Templates/OPLS2005/HeteroAtoms/" ]; then
		mkdir -p "${Pele_directory}/DataLocal/Templates/OPLS2005/HeteroAtoms/"
	fi
	if [ ! -d "${Pele_directory}/DataLocal/Templates/AMBER99sbBSC0/HeteroAtoms/" ]; then
		mkdir -p "${Pele_directory}/DataLocal/Templates/AMBER99sbBSC0/HeteroAtoms/"
	fi
	if [ ! -d "${Pele_directory}/DataLocal/Templates/AMBER99/HeteroAtoms/" ]; then
		mkdir -p "${Pele_directory}/DataLocal/Templates/AMBER99sb/HeteroAtoms/"
	fi

	if [[ "$FORCEFIELD" == "AMBER99sbBSC0" ]]; then
		mv "${output_files[0]}" "${Pele_directory}/DataLocal/Templates/AMBER99sbBSC0/HeteroAtoms/"
	elif [[ "$FORCEFIELD" == "AMBER99sb" ]]; then
		mv "${output_files[0]}" "${Pele_directory}/DataLocal/Templates/AMBER99sb/HeteroAtoms/"
	else
		mv "${output_files[0]}" "${Pele_directory}/DataLocal/Templates/OPLS2005/HeteroAtoms/"
	fi
	
	
	if [ ! -d "${Pele_directory}/DataLocal/LigandRotamerLibs" ]; then
		mkdir -p "${Pele_directory}/DataLocal/LigandRotamerLibs"
	fi
	mv "${output_files[1]}" "${Pele_directory}/DataLocal/LigandRotamerLibs"


	if [ ! -d "${Pele_directory}/results" ]; then
		mkdir "${Pele_directory}/results"
	else
		rm -rf "${Pele_directory}/results"
		mkdir "${Pele_directory}/results"
	fi

	cp $PDB_FILE "${Pele_directory}/complex.pdb"
	cp $NATIVE "${Pele_directory}/native.pdb"

	if [[ -f "$CONFILE" ]]; then
		cp "$CONFILE" "${Pele_directory}/control_file"
	else
		cp ${PlopRotTemp}/PeleTemplates/control_file "${Pele_directory}/control_file"
	fi
	
	sed -i 's,$CHAIN,'"${LIGAND_CHAIN}"',g' "${Pele_directory}/control_file"
	sed -i 's,$FORCEFIELD,'"${FORCEFIELD}"',g' "${Pele_directory}/control_file"

	#sed -i 's,$CHAIN,'"${LIGAND_CHAIN}"',g' "${Pele_directory}/pca_control_file"

	#rm ligand_mae
	rm input.txt

	#prepare input
	python Helpers/input_prep.py ${Pele_directory}/complex.pdb

	echo "input prepared"

	cd "$Pele_directory"
	ln -s /home/dani/repos/PELErev12535/Data/ Data
	ln -s /home/dani/repos/PELErev12535/Documents/ Documents

	#RunPele
	mpirun -np $CPUS PELE-1.5_mpi control_file --license-directory ~/repos/PELErev12535/licenses/

fi

###################################################################################







