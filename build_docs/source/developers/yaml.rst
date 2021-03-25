YAML Model
==========


*******************************************************
General settings
*******************************************************


Package
########################################################################






:YAML key: package


:Parser key: package
:Type: ``str``






Processors:
 * ``str_validator`` from `Pydantic`




Adaptive
########################################################################






:YAML key: adaptive


:Parser key: adaptive
:Type: ``str``






Processors:
 * ``str_validator`` from `Pydantic`

 * ``set_adaptive`` from `SimulationParamsModel`




Global Inputs / Input
########################################################################






:YAML key: global_inputs


:Parser key: input
:Type: ``Any``






Processors:
 * ``set_input_glob`` from `SimulationParamsModel`




System
########################################################################






:YAML key: system


:Parser key: system
:Type: ``str``

:Default: ``''``






Processors:
 * ``str_validator`` from `Pydantic`

 * ``construct_path`` from `YamlParserModel`

 * ``set_system_glob`` from `SimulationParamsModel`

 * ``validate_adaptive_required_fields`` from `SimulationParamsModel`




Resname / Residue
########################################################################






:YAML key: resname


:Parser key: residue
:Type: ``str``






Processors:
 * ``str_validator`` from `Pydantic`

 * ``validate_residue`` from `YamlParserModel`

 * ``validate_adaptive_required_fields`` from `SimulationParamsModel`




Chain
########################################################################






:YAML key: chain


:Parser key: chain
:Type: ``str``






Processors:
 * ``str_validator`` from `Pydantic`

 * ``validate_adaptive_required_fields`` from `SimulationParamsModel`




Test
########################################################################






:YAML key: test


:Parser key: test
:Type: ``bool``






Processors:
 * ``bool_validator`` from `Pydantic`




Verbose
########################################################################






:YAML key: verbose


:Parser key: verbose
:Type: ``bool``






Processors:
 * ``bool_validator`` from `Pydantic`

 * ``format_verbose`` from `SimulationParamsModel`




Debug
########################################################################






:YAML key: debug


:Parser key: debug
:Type: ``bool``

:Default: ``False``






Processors:
 * ``bool_validator`` from `Pydantic`




Working Folder / Folder
########################################################################






:YAML key: working_folder


:Parser key: folder
:Type: ``str``






Processors:
 * ``str_validator`` from `Pydantic`




Output
########################################################################






:YAML key: output


:Parser key: output
:Type: ``str``

:Default: ``'output'``






Processors:
 * ``str_validator`` from `Pydantic`





*******************************************************
Other
*******************************************************


Hbond
########################################################################




.. warning::
    This is a candidate for deprecation.



:YAML key: hbond


:Parser key: hbond
:Type: ``Any``

:Default: ``[None, None]``








Pele
########################################################################




.. warning::
    This is a candidate for deprecation.



:YAML key: pele


:Parser key: pele
:Type: ``Any``








Temp
########################################################################




.. warning::
    This is a candidate for deprecation.



:YAML key: None, value from `temperature` parser


:Parser key: temp
:Type: ``int``






Processors:
 * ``int_validator`` from `Pydantic`




Adaptive Restart
########################################################################






:YAML key: adaptive_restart


:Parser key: adaptive_restart
:Type: ``bool``






Processors:
 * ``bool_validator`` from `Pydantic`




Templates / Template
########################################################################






:YAML key: templates


:Parser key: template
:Type: ``str``








Ext Rotamers
########################################################################






:YAML key: None, value from `rotamers` parser


:Parser key: ext_rotamers
:Type: ``str``








TERs / Gaps Ter
########################################################################






:YAML key: TERs


:Parser key: gaps_ter
:Type: ``bool``




If value is ``None``,
it will fall back to ``gaps_ter`` software setting.
If there is no software setting,
it will default to ``False``.




Processors:
 * ``bool_validator`` from `Pydantic`




Charge Ters / Charge Ter
########################################################################






:YAML key: charge_ters


:Parser key: charge_ter
:Type: ``bool``




If value is ``None``,
it will fall back to ``charge_ter`` software setting.
If there is no software setting,
it will default to ``False``.




Processors:
 * ``bool_validator`` from `Pydantic`




Mpi Params
########################################################################






:YAML key: mpi_params


:Parser key: mpi_params
:Type: ``Any``






Processors:
 * ``format_mpi_params`` from `SimulationParamsModel`




Nonstandard
########################################################################






:YAML key: nonstandard


:Parser key: nonstandard
:Type: ``str``




If value is  falsy,
it will fall back to ``nonstandard`` software setting.
If there is no software setting,
it will default to ``[]``.






Prepwizard
########################################################################






:YAML key: prepwizard


:Parser key: prepwizard
:Type: ``bool``






Processors:
 * ``bool_validator`` from `Pydantic`




proximityDetection
########################################################################






:YAML key: proximityDetection


:Parser key: proximityDetection
:Type: ``bool``






Processors:
 * ``bool_validator`` from `Pydantic`

 * ``set_proximityDetection`` from `SimulationParamsModel`




Poses
########################################################################






:YAML key: poses


:Parser key: poses
:Type: ``int``






Processors:
 * ``int_validator`` from `Pydantic`




Precision Glide
########################################################################




.. warning::
    This is a candidate for deprecation.



:YAML key: precision_glide


:Parser key: precision_glide
:Type: ``Any``








Msm
########################################################################






:YAML key: msm


:Parser key: msm
:Type: ``Any``








Precision
########################################################################




.. warning::
    This is a candidate for deprecation.



:YAML key: precision


:Parser key: precision
:Type: ``Any``








Exit Clust / Clust
########################################################################




.. warning::
    This is a candidate for deprecation.



:YAML key: exit_clust


:Parser key: clust
:Type: ``Any``


:Tests value: ``2``







Restart
########################################################################






:YAML key: restart


:Parser key: restart
:Type: ``Any``




If value is  falsy,
it will fall back to ``restart`` software setting.
If there is no software setting,
it will default to ``'all'``.






Lagtime
########################################################################




.. warning::
    This is a candidate for deprecation.



:YAML key: lagtime


:Parser key: lagtime
:Type: ``Any``








Msm Clust
########################################################################






:YAML key: msm_clust


:Parser key: msm_clust
:Type: ``Any``








Rescoring
########################################################################






:YAML key: rescoring


:Parser key: rescoring
:Type: ``bool``






Processors:
 * ``bool_validator`` from `Pydantic`




In Out
########################################################################






:YAML key: in_out


:Parser key: in_out
:Type: ``bool``






Processors:
 * ``bool_validator`` from `Pydantic`




In Out Soft
########################################################################






:YAML key: in_out_soft


:Parser key: in_out_soft
:Type: ``bool``






Processors:
 * ``bool_validator`` from `Pydantic`




Exit
########################################################################






:YAML key: exit


:Parser key: exit
:Type: ``Any``








Exit Value
########################################################################






:YAML key: exit_value


:Parser key: exit_value
:Type: ``float``






Processors:
 * ``float_validator`` from `Pydantic`




Exit Condition
########################################################################






:YAML key: exit_condition


:Parser key: exit_condition
:Type: ``str``






Processors:
 * ``str_validator`` from `Pydantic`




Exit Trajnum
########################################################################






:YAML key: exit_trajnum


:Parser key: exit_trajnum
:Type: ``int``






Processors:
 * ``int_validator`` from `Pydantic`




One Exit
########################################################################






:YAML key: one_exit


:Parser key: one_exit
:Type: ``Any``








Box Type
########################################################################






:YAML key: box_type


:Parser key: box_type
:Type: ``Any``








Box Metric
########################################################################






:YAML key: box_metric


:Parser key: box_metric
:Type: ``Any``








Time
########################################################################






:YAML key: time


:Parser key: time
:Type: ``Any``








Nosasa
########################################################################






:YAML key: nosasa


:Parser key: nosasa
:Type: ``Any``








Perc Sasa
########################################################################






:YAML key: perc_sasa


:Parser key: perc_sasa
:Type: ``Any``








Seed
########################################################################






:YAML key: seed


:Parser key: seed
:Type: ``int``

:Default: ``generate_random_seed()``






Processors:
 * ``int_validator`` from `Pydantic`




Pdb
########################################################################






:YAML key: pdb


:Parser key: pdb
:Type: ``bool``






Processors:
 * ``bool_validator`` from `Pydantic`

 * ``check_extensions`` from `SimulationParamsModel`




Log
########################################################################






:YAML key: log


:Parser key: log
:Type: ``Any``








Nonrenum
########################################################################






:YAML key: nonrenum


:Parser key: nonrenum
:Type: ``Any``








Pele Exec
########################################################################






:YAML key: pele_exec


:Parser key: pele_exec
:Type: ``str``

:Default: ``'bin/Pele_mpi'``






Processors:
 * ``str_validator`` from `Pydantic`




Pele Data
########################################################################






:YAML key: pele_data


:Parser key: pele_data
:Type: ``str``

:Default: ``'Data'``






Processors:
 * ``str_validator`` from `Pydantic`




Pele Documents
########################################################################






:YAML key: pele_documents


:Parser key: pele_documents
:Type: ``str``

:Default: ``'Documents'``






Processors:
 * ``str_validator`` from `Pydantic`




Pca
########################################################################






:YAML key: pca


:Parser key: pca
:Type: ``Any``






Processors:
 * ``format_pca`` from `SimulationParamsModel`




Anm Direction
########################################################################






:YAML key: anm_direction


:Parser key: anm_direction
:Type: ``str``




If value is  falsy,
it will fall back to ``anm_direction`` software setting.
If there is no software setting,
it will default to ``'random'``.




Processors:
 * ``str_validator`` from `Pydantic`




Anm Mix Modes
########################################################################






:YAML key: anm_mix_modes


:Parser key: anm_mix_modes
:Type: ``str``




If value is  falsy,
it will fall back to ``anm_mix_modes`` software setting.
If there is no software setting,
it will default to ``'mixMainModeWithOthersModes'``.




Processors:
 * ``str_validator`` from `Pydantic`




Anm Picking Mode
########################################################################






:YAML key: anm_picking_mode


:Parser key: anm_picking_mode
:Type: ``str``




If value is  falsy,
it will fall back to ``anm_picking_mode`` software setting.
If there is no software setting,
it will default to ``'RANDOM_MODE'``.




Processors:
 * ``str_validator`` from `Pydantic`




Anm Displacement
########################################################################






:YAML key: anm_displacement


:Parser key: anm_displacement
:Type: ``float``




If value is  falsy,
it will fall back to ``anm_displacement`` software setting.
If there is no software setting,
it will default to ``0.75``.




Processors:
 * ``float_validator`` from `Pydantic`




Anm Modes Change
########################################################################






:YAML key: anm_modes_change


:Parser key: anm_modes_change
:Type: ``Any``




If value is  falsy,
it will fall back to ``anm_modes_change`` software setting.
If there is no software setting,
it will default to ``4``.






Anm Num Of Modes
########################################################################






:YAML key: anm_num_of_modes


:Parser key: anm_num_of_modes
:Type: ``int``




If value is  falsy,
it will fall back to ``anm_num_of_modes`` software setting.
If there is no software setting,
it will default to ``6``.




Processors:
 * ``int_validator`` from `Pydantic`




Anm Relaxation Constr
########################################################################






:YAML key: anm_relaxation_constr


:Parser key: anm_relaxation_constr
:Type: ``float``




If value is  falsy,
it will fall back to ``anm_relaxation_constr`` software setting.
If there is no software setting,
it will default to ``0.5``.




Processors:
 * ``float_validator`` from `Pydantic`




Remove Constraints
########################################################################






:YAML key: remove_constraints


:Parser key: remove_constraints
:Type: ``bool``




If value is ``None``,
it will fall back to ``remove_constraints`` software setting.
If there is no software setting,
it will default to ``False``.




Processors:
 * ``bool_validator`` from `Pydantic`




Pca Traj
########################################################################






:YAML key: pca_traj


:Parser key: pca_traj
:Type: ``str``








Perturbation
########################################################################






:YAML key: perturbation


:Parser key: perturbation
:Type: ``Any``






Processors:
 * ``set_perturbation`` from `SimulationParamsModel`




Sasa
########################################################################






:YAML key: sasa


:Parser key: sasa
:Type: ``str``




If value is  falsy,
it will fall back to ``sasa`` software setting.
If there is no software setting,
it will default to ``'\n                        { "type": "sasa",\n\n                           "tag": "sasaLig",\n\n                           "selection": { "chains": { "names": ["$CHAIN"] } }\n\n                        },\n'``.




Processors:
 * ``str_validator`` from `Pydantic`

 * ``only_with_perturbation`` from `SimulationParamsModel`




Binding Energy
########################################################################






:YAML key: binding_energy


:Parser key: binding_energy
:Type: ``str``




If value is  falsy,
it will fall back to ``binding_energy`` software setting.
If there is no software setting,
it will default to ``'\n                        { "type": "bindingEnergy",\n\n                           "boundPartSelection": { "chains": { "names": ["$CHAIN"] } }\n\n                        },\n'``.




Processors:
 * ``str_validator`` from `Pydantic`

 * ``only_with_perturbation`` from `SimulationParamsModel`




Parameters
########################################################################






:YAML key: parameters


:Parser key: parameters
:Type: ``str``




If value is  falsy,
it will fall back to ``params`` software setting.
If there is no software setting,
it will default to ``True``.




Processors:
 * ``str_validator`` from `Pydantic`

 * ``only_with_perturbation`` from `SimulationParamsModel`




Analyse
########################################################################






:YAML key: analyse


:Parser key: analyse
:Type: ``Any``

:Default: ``True``








Selection To Perturb
########################################################################






:YAML key: selection_to_perturb


:Parser key: selection_to_perturb
:Type: ``str``




If value is  falsy,
it will fall back to ``selection_to_perturb`` software setting.
If there is no software setting,
it will default to ``'"selectionToPerturb" : { "chains" : { "names" : [ "$CHAIN" ] } },'``.




Processors:
 * ``str_validator`` from `Pydantic`

 * ``only_with_perturbation`` from `SimulationParamsModel`




Mae
########################################################################






:YAML key: mae


:Parser key: mae
:Type: ``Any``

:Default: ``False``








Constrain Core
########################################################################






:YAML key: constrain_core


:Parser key: constrain_core
:Type: ``Any``




If value is  falsy,
it will fall back to ``constrain_core`` software setting.






Spawning Condition
########################################################################






:YAML key: spawning_condition


:Parser key: spawning_condition
:Type: ``Any``




If value is  falsy,
it will fall back to ``spawning_condition`` software setting.




Processors:
 * ``format_spawning_condition`` from `SimulationParamsModel`




Overwrite Analysis / Overwrite
########################################################################






:YAML key: overwrite_analysis


:Parser key: overwrite
:Type: ``bool``

:Default: ``True``






Processors:
 * ``bool_validator`` from `Pydantic`




Analysis Nclust
########################################################################






:YAML key: analysis_nclust


:Parser key: analysis_nclust
:Type: ``int``

:Default: ``10``






Processors:
 * ``int_validator`` from `Pydantic`




Te Column
########################################################################






:YAML key: te_column


:Parser key: te_column
:Type: ``int``

:Default: ``4``






Processors:
 * ``int_validator`` from `Pydantic`




Be Column
########################################################################






:YAML key: be_column


:Parser key: be_column
:Type: ``int``

:Default: ``5``






Processors:
 * ``int_validator`` from `Pydantic`




Limit Column
########################################################################






:YAML key: limit_column


:Parser key: limit_column
:Type: ``int``

:Default: ``6``






Processors:
 * ``int_validator`` from `Pydantic`




COMligandConstraint / Com
########################################################################






:YAML key: COMligandConstraint


:Parser key: com
:Type: ``Any``




If value is  falsy,
it will fall back to ``COMligandConstraint`` software setting.
If there is no software setting,
it will default to ``0``.






Pele License
########################################################################






:YAML key: pele_license


:Parser key: pele_license
:Type: ``str``

:Default: ``'licenses'``






Processors:
 * ``str_validator`` from `Pydantic`




License
########################################################################






:YAML key: None, value from `pele_license` parser


:Parser key: license
:Type: ``Any``








Schrodinger
########################################################################






:YAML key: schrodinger


:Parser key: schrodinger
:Type: ``str``






Processors:
 * ``str_validator`` from `Pydantic`




No Check
########################################################################






:YAML key: no_check


:Parser key: no_check
:Type: ``bool``

:Default: ``False``






Processors:
 * ``bool_validator`` from `Pydantic`




Frag Pele
########################################################################





:YAML key: None, value calculated in simulation params

:Parser key: frag_pele
:Type: ``Any``






Processors:
 * ``set_frag_pele`` from `SimulationParamsModel`




Complexes
########################################################################





:YAML key: None, value calculated in simulation params

:Parser key: complexes
:Type: ``Any``






Processors:
 * ``set_complexes`` from `SimulationParamsModel`




Frag Pele Steps
########################################################################





:YAML key: None, value calculated in simulation params

:Parser key: frag_pele_steps
:Type: ``Any``






Processors:
 * ``set_frag_pele_steps`` from `SimulationParamsModel`




Output Path
########################################################################





:YAML key: None, value calculated in simulation params

:Parser key: output_path
:Type: ``Any``






Processors:
 * ``set_output_path`` from `SimulationParamsModel`




Logfile
########################################################################





:YAML key: None, value calculated in simulation params

:Parser key: logfile
:Type: ``str``






Processors:
 * ``str_validator`` from `Pydantic`

 * ``set_logfile`` from `SimulationParamsModel`




Water
########################################################################





:YAML key: None, value calculated in simulation params

:Parser key: water
:Type: ``str``

:Default: ``''``






Processors:
 * ``str_validator`` from `Pydantic`




Ligand
########################################################################





:YAML key: None, value calculated in simulation params

:Parser key: ligand
:Type: ``Any``

:Default: ``'"ligandResname" : "$LIG_RES",'``






Processors:
 * ``only_with_perturbation`` from `SimulationParamsModel`




External Template
########################################################################





:YAML key: None, value calculated in simulation params

:Parser key: external_template
:Type: ``Any``




If value is  falsy,
it will fall back to ``template`` software setting.
If there is no software setting,
it will default to ``[]``.






External Rotamers
########################################################################





:YAML key: None, value calculated in simulation params

:Parser key: external_rotamers
:Type: ``Any``




If value is  falsy,
it will fall back to ``rotamers`` software setting.
If there is no software setting,
it will default to ``[]``.






Spython
########################################################################





:YAML key: None, value calculated in simulation params

:Parser key: spython
:Type: ``Any``

:Default: ``'utilities/python'``






Processors:
 * ``check_spython_path`` from `SimulationParamsModel`




Lig
########################################################################





:YAML key: None, value calculated in simulation params

:Parser key: lig
:Type: ``Any``






Processors:
 * ``set_lig`` from `SimulationParamsModel`




Sasa Max
########################################################################





:YAML key: None, value calculated in simulation params

:Parser key: sasa_max
:Type: ``Any``








Sasa Min
########################################################################





:YAML key: None, value calculated in simulation params

:Parser key: sasa_min
:Type: ``Any``








Clust / Clusters
########################################################################





:YAML key: None, value calculated in simulation params

:Parser key: clusters
:Type: ``Any``


:Tests value: ``2``







Allow Empty Selectors
########################################################################





:YAML key: None, value calculated in simulation params

:Parser key: allow_empty_selectors
:Type: ``Any``






Processors:
 * ``format_allow_empty_selectors`` from `SimulationParamsModel`




Templates
########################################################################





:YAML key: None, value calculated in simulation params

:Parser key: templates
:Type: ``Any``

:Default: ``'/Users/agruszka/Projects/pele_platform/pele_platform/PeleTemplates'``








Xtc
########################################################################





:YAML key: None, value calculated in simulation params

:Parser key: xtc
:Type: ``bool``






Processors:
 * ``bool_validator`` from `Pydantic`

 * ``check_extensions`` from `SimulationParamsModel`




Constraints
########################################################################





:YAML key: None, value calculated in simulation params

:Parser key: constraints
:Type: ``Any``








Water Energy
########################################################################





:YAML key: None, value calculated in simulation params

:Parser key: water_energy
:Type: ``Any``









*******************************************************
Simulation parameters
*******************************************************


Forcefield
########################################################################






:YAML key: forcefield


:Parser key: forcefield
:Type: ``str``

:Default: ``'OPLS2005'``




If value is  falsy,
it will fall back to ``forcefield`` software setting.
If there is no software setting,
it will default to ``'OPLS2005'``.




Processors:
 * ``str_validator`` from `Pydantic`




Anm Freq
########################################################################






:YAML key: anm_freq


:Parser key: anm_freq
:Type: ``int``


:Tests value: ``0``



If value is ``None``,
it will fall back to ``anm_freq`` software setting.
If there is no software setting,
it will default to ``4``.




Processors:
 * ``int_validator`` from `Pydantic`




Sidechain Freq
########################################################################






:YAML key: sidechain_freq


:Parser key: sidechain_freq
:Type: ``int``


:Tests value: ``0``



If value is ``None``,
it will fall back to ``sidechain_freq`` software setting.
If there is no software setting,
it will default to ``2``.




Processors:
 * ``int_validator`` from `Pydantic`




Min Freq
########################################################################






:YAML key: min_freq


:Parser key: min_freq
:Type: ``int``


:Tests value: ``0``



If value is ``None``,
it will fall back to ``min_freq`` software setting.
If there is no software setting,
it will default to ``1``.




Processors:
 * ``int_validator`` from `Pydantic`




Water Freq
########################################################################






:YAML key: water_freq


:Parser key: water_freq
:Type: ``int``




If value is ``None``,
it will fall back to ``water_freq`` software setting.
If there is no software setting,
it will default to ``1``.




Processors:
 * ``int_validator`` from `Pydantic`




Temperature
########################################################################






:YAML key: temperature


:Parser key: temperature
:Type: ``int``


:Tests value: ``10000``



If value is  falsy,
it will fall back to ``temperature`` software setting.
If there is no software setting,
it will default to ``1500``.




Processors:
 * ``int_validator`` from `Pydantic`




Sidechain Res / Sidechain Resolution
########################################################################






:YAML key: sidechain_res


:Parser key: sidechain_resolution
:Type: ``int``




If value is  falsy,
it will fall back to ``sidechain_resolution`` software setting.
If there is no software setting,
it will default to ``30``.




Processors:
 * ``int_validator`` from `Pydantic`




Steric Trials
########################################################################






:YAML key: steric_trials


:Parser key: steric_trials
:Type: ``int``




If value is  falsy,
it will fall back to ``steric_trials`` software setting.
If there is no software setting,
it will default to ``250``.




Processors:
 * ``int_validator`` from `Pydantic`




Overlap Factor
########################################################################






:YAML key: overlap_factor


:Parser key: overlap_factor
:Type: ``float``




If value is  falsy,
it will fall back to ``overlap_factor`` software setting.
If there is no software setting,
it will default to ``0.65``.




Processors:
 * ``float_validator`` from `Pydantic`




Steering
########################################################################

Number of translations in the same direction.




:YAML key: steering


:Parser key: steering
:Type: ``int``




If value is  falsy,
it will fall back to ``steering`` software setting.
If there is no software setting,
it will default to ``0``.




Processors:
 * ``int_validator`` from `Pydantic`




Solvent
########################################################################






:YAML key: solvent


:Parser key: solvent
:Type: ``str``




If value is  falsy,
it will fall back to ``solvent`` software setting.
If there is no software setting,
it will default to ``'VDGBNP'``.




Processors:
 * ``str_validator`` from `Pydantic`




Spawning
########################################################################






:YAML key: spawning


:Parser key: spawning
:Type: ``str``




If value is  falsy,
it will fall back to ``spawning_type`` software setting.
If there is no software setting,
it will default to ``'independent'``.




Processors:
 * ``str_validator`` from `Pydantic`




Iterations
########################################################################






:YAML key: iterations


:Parser key: iterations
:Type: ``int``


:Tests value: ``1``



If value is  falsy,
it will fall back to ``iterations`` software setting.




Processors:
 * ``int_validator`` from `Pydantic`

 * ``set_iterations`` from `SimulationParamsModel`




Steps
########################################################################




.. warning::
    This is a candidate for deprecation.



:YAML key: steps


:Parser key: steps
:Type: ``int``


:Tests value: ``1``





Processors:
 * ``int_validator`` from `Pydantic`




Pele Steps
########################################################################






:YAML key: None, value from `steps` parser


:Parser key: pele_steps
:Type: ``int``




If value is  falsy,
it will fall back to ``pele_steps`` software setting.
If there is no software setting,
it will default to ``8``.




Processors:
 * ``int_validator`` from `Pydantic`




Cpus
########################################################################






:YAML key: cpus


:Parser key: cpus
:Type: ``int``


:Tests value: ``5``



If value is  falsy,
it will fall back to ``cpus`` software setting.
If there is no software setting,
it will default to ``60``.




Processors:
 * ``int_validator`` from `Pydantic`




Density
########################################################################






:YAML key: density


:Parser key: density
:Type: ``str``




If value is  falsy,
it will fall back to ``density`` software setting.
If there is no software setting,
it will default to ``'null'``.




Processors:
 * ``str_validator`` from `Pydantic`




Cluster Values
########################################################################






:YAML key: cluster_values


:Parser key: cluster_values
:Type: ``float``




If value is  falsy,
it will fall back to ``cluster_values`` software setting.
If there is no software setting,
it will default to ``'[1.75, 2.5, 4, 6]'``.






Cluster Conditions
########################################################################






:YAML key: cluster_conditions


:Parser key: cluster_conditions
:Type: ``float``




If value is  falsy,
it will fall back to ``cluster_conditions`` software setting.
If there is no software setting,
it will default to ``'[1, 0.6, 0.4, 0.0]'``.






Simulation Type
########################################################################






:YAML key: simulation_type


:Parser key: simulation_type
:Type: ``str``




If value is  falsy,
it will fall back to ``simulation_type`` software setting.
If there is no software setting,
it will default to ``'pele'``.




Processors:
 * ``str_validator`` from `Pydantic`




Equilibration
########################################################################






:YAML key: equilibration


:Parser key: equilibration
:Type: ``bool``






Processors:
 * ``bool_validator`` from `Pydantic`

 * ``format_verbose`` from `SimulationParamsModel`




Clust Type
########################################################################






:YAML key: clust_type


:Parser key: clust_type
:Type: ``str``




If value is  falsy,
it will fall back to ``clust_type`` software setting.
If there is no software setting,
it will default to ``'rmsd'``.




Processors:
 * ``str_validator`` from `Pydantic`




Equilibration Steps / Eq Steps
########################################################################






:YAML key: equilibration_steps


:Parser key: eq_steps
:Type: ``int``






Processors:
 * ``int_validator`` from `Pydantic`




Report / Report Name
########################################################################




.. warning::
    This is a candidate for deprecation.



:YAML key: report


:Parser key: report_name
:Type: ``str``

:Default: ``'report'``






Processors:
 * ``str_validator`` from `Pydantic`




Traj / Traj Name
########################################################################






:YAML key: traj


:Parser key: traj_name
:Type: ``str``

:Default: ``'trajectory.pdb'``






Processors:
 * ``str_validator`` from `Pydantic`




Epsilon
########################################################################






:YAML key: epsilon


:Parser key: epsilon
:Type: ``float``




If value is  falsy,
it will fall back to ``epsilon`` software setting.
If there is no software setting,
it will default to ``0``.




Processors:
 * ``float_validator`` from `Pydantic`




Bias Column
########################################################################






:YAML key: bias_column


:Parser key: bias_column
:Type: ``int``




If value is  falsy,
it will fall back to ``bias_column`` software setting.
If there is no software setting,
it will default to ``5``.




Processors:
 * ``int_validator`` from `Pydantic`




Randomize
########################################################################






:YAML key: randomize


:Parser key: randomize
:Type: ``bool``




If value is  falsy,
it will fall back to ``randomize`` software setting.
If there is no software setting,
it will default to ``False``.




Processors:
 * ``bool_validator`` from `Pydantic`





*******************************************************
Miscellaneous
*******************************************************


Usesrun
########################################################################






:YAML key: usesrun


:Parser key: usesrun
:Type: ``bool``






Processors:
 * ``bool_validator`` from `Pydantic`

 * ``set_usesrun`` from `SimulationParamsModel`





*******************************************************
Out in
*******************************************************


Out In
########################################################################






:YAML key: out_in


:Parser key: out_in
:Type: ``bool``






Processors:
 * ``bool_validator`` from `Pydantic`




Initial Site
########################################################################






:YAML key: initial_site


:Parser key: initial_site
:Type: ``str``






Processors:
 * ``str_validator`` from `Pydantic`

 * ``validate_atom_string`` from `YamlParserModel`





*******************************************************
Ligand preparation
*******************************************************


Gridres
########################################################################






:YAML key: gridres


:Parser key: gridres
:Type: ``int``

:Default: ``10``






Processors:
 * ``int_validator`` from `Pydantic`




Core
########################################################################






:YAML key: core


:Parser key: core
:Type: ``int``

:Default: ``-1``






Processors:
 * ``int_validator`` from `Pydantic`




Maxtorsion / Mtor
########################################################################






:YAML key: maxtorsion


:Parser key: mtor
:Type: ``int``

:Default: ``4``






Processors:
 * ``int_validator`` from `Pydantic`




N
########################################################################

Maximum number of flexible sidechains in the ligand.




:YAML key: n


:Parser key: n
:Type: ``int``

:Default: ``10000``






Processors:
 * ``int_validator`` from `Pydantic`




Ext Temp
########################################################################






:YAML key: None, value from `template` parser


:Parser key: ext_temp
:Type: ``str``








Rotamers
########################################################################






:YAML key: rotamers


:Parser key: rotamers
:Type: ``str``








Mae Lig
########################################################################

Maestro file to extract quantum charges.




:YAML key: mae_lig


:Parser key: mae_lig
:Type: ``str``




If value is  falsy,
it will fall back to ``mae_lig`` software setting.




Processors:
 * ``str_validator`` from `Pydantic`

 * ``construct_path`` from `YamlParserModel`




Skip Ligand Prep
########################################################################






:YAML key: skip_ligand_prep


:Parser key: skip_ligand_prep
:Type: ``str``




If value is  falsy,
it will fall back to ``args.skip_ligand_prep`` software setting.
If there is no software setting,
it will default to ``[]``.







*******************************************************
Protein preparation
*******************************************************


Skip Preprocess / No Ppp
########################################################################






:YAML key: skip_preprocess


:Parser key: no_ppp
:Type: ``bool``




If value is  falsy,
it will fall back to ``no_ppp`` software setting.
If there is no software setting,
it will default to ``False``.




Processors:
 * ``bool_validator`` from `Pydantic`




Skip Preprocess / Skip Prep
########################################################################






:YAML key: skip_preprocess


:Parser key: skip_prep
:Type: ``bool``




If value is  falsy,
it will fall back to ``skip_prep`` software setting.
If there is no software setting,
it will default to ``False``.




Processors:
 * ``bool_validator`` from `Pydantic`





*******************************************************
Box settings
*******************************************************


Box Radius
########################################################################






:YAML key: box_radius


:Parser key: box_radius
:Type: ``float``




If value is  falsy,
it will fall back to ``box_radius`` software setting.




Processors:
 * ``float_validator`` from `Pydantic`




Box Center
########################################################################






:YAML key: box_center


:Parser key: box_center
:Type: ``Union[List[float], str]``






Processors:
 * ``calculate_box_center`` from `SimulationParamsModel`




Box
########################################################################






:YAML key: box


:Parser key: box
:Type: ``Any``









*******************************************************
Metrics
*******************************************************


Rmsd Pdb / Native
########################################################################






:YAML key: rmsd_pdb


:Parser key: native
:Type: ``str``

:Default: ``False``






Processors:
 * ``str_validator`` from `Pydantic`




Atom Dist
########################################################################






:YAML key: atom_dist


:Parser key: atom_dist
:Type: ``Union[List[str], List[int]]``

:Default: ``list()``









*******************************************************
Global Exploration
*******************************************************


Randomize
########################################################################






:YAML key: randomize


:Parser key: randomize
:Type: ``bool``




If value is  falsy,
it will fall back to ``randomize`` software setting.
If there is no software setting,
it will default to ``False``.




Processors:
 * ``bool_validator`` from `Pydantic`




Global / Full
########################################################################






:YAML key: global


:Parser key: full
:Type: ``bool``






Processors:
 * ``bool_validator`` from `Pydantic`





*******************************************************
Water
*******************************************************


Waters
########################################################################






:YAML key: waters


:Parser key: waters
:Type: ``Union[str, List[str]]``




If value is  falsy,
it will fall back to ``waters`` software setting.
If there is no software setting,
it will default to ``[]``.






Water Center
########################################################################






:YAML key: water_center


:Parser key: water_center
:Type: ``Union[List[float], str]``








Water Temp
########################################################################






:YAML key: water_temp


:Parser key: water_temp
:Type: ``Any``




If value is  falsy,
it will fall back to ``water_temp`` software setting.
If there is no software setting,
it will default to ``5000``.






Water Overlap
########################################################################






:YAML key: water_overlap


:Parser key: water_overlap
:Type: ``float``




If value is  falsy,
it will fall back to ``water_overlap`` software setting.
If there is no software setting,
it will default to ``0.78``.




Processors:
 * ``float_validator`` from `Pydantic`




Water Constr
########################################################################






:YAML key: water_constr


:Parser key: water_constr
:Type: ``Any``




If value is  falsy,
it will fall back to ``water_constr`` software setting.
If there is no software setting,
it will default to ``0``.






Water Trials
########################################################################






:YAML key: water_trials


:Parser key: water_trials
:Type: ``Any``




If value is  falsy,
it will fall back to ``water_trials`` software setting.
If there is no software setting,
it will default to ``10000``.






Water Radius
########################################################################






:YAML key: water_radius


:Parser key: water_radius
:Type: ``int``

:Default: ``6``






Processors:
 * ``int_validator`` from `Pydantic`




Water Empty Selector
########################################################################






:YAML key: water_empty_selector


:Parser key: water_empty_selector
:Type: ``Any``

:Default: ``False``








N Waters
########################################################################






:YAML key: n_waters


:Parser key: n_waters
:Type: ``int``




If value is  falsy,
it will fall back to ``n_waters`` software setting.
If there is no software setting,
it will default to ``0``.




Processors:
 * ``int_validator`` from `Pydantic`





*******************************************************
Induced fit
*******************************************************


Induced Fit Exhaustive
########################################################################






:YAML key: induced_fit_exhaustive


:Parser key: induced_fit_exhaustive
:Type: ``bool``






Processors:
 * ``bool_validator`` from `Pydantic`




Induced Fit Fast
########################################################################






:YAML key: induced_fit_fast


:Parser key: induced_fit_fast
:Type: ``bool``






Processors:
 * ``bool_validator`` from `Pydantic`





*******************************************************
FragPELE
*******************************************************


Frag
########################################################################






:YAML key: frag


:Parser key: frag
:Type: ``Any``








Cleanup
########################################################################

Automatically cleans up fragment files, only applicable to FragPELE.




:YAML key: cleanup


:Parser key: cleanup
:Type: ``bool``

:Default: ``False``






Processors:
 * ``bool_validator`` from `Pydantic`




Frag Run
########################################################################






:YAML key: frag_run


:Parser key: frag_run
:Type: ``bool``

:Default: ``True``






Processors:
 * ``bool_validator`` from `Pydantic`




Frag Core
########################################################################






:YAML key: frag_core


:Parser key: frag_core
:Type: ``str``






Processors:
 * ``str_validator`` from `Pydantic`




Frag Input
########################################################################






:YAML key: frag_input


:Parser key: frag_input
:Type: ``str``

:Default: ``False``






Processors:
 * ``str_validator`` from `Pydantic`




Frag Ligands
########################################################################






:YAML key: frag_ligands


:Parser key: frag_ligands
:Type: ``Any``

:Default: ``False``








Growing Steps
########################################################################






:YAML key: growing_steps


:Parser key: growing_steps
:Type: ``int``

:Default: ``False``






Processors:
 * ``int_validator`` from `Pydantic`




Steps In Gs / Frag Steps
########################################################################






:YAML key: steps_in_gs


:Parser key: frag_steps
:Type: ``int``

:Default: ``False``






Processors:
 * ``int_validator`` from `Pydantic`




Sampling Steps / Frag Eq Steps
########################################################################






:YAML key: sampling_steps


:Parser key: frag_eq_steps
:Type: ``int``

:Default: ``False``






Processors:
 * ``int_validator`` from `Pydantic`




Protocol
########################################################################






:YAML key: protocol


:Parser key: protocol
:Type: ``str``






Processors:
 * ``str_validator`` from `Pydantic`




Frag Ai
########################################################################






:YAML key: frag_ai


:Parser key: frag_ai
:Type: ``bool``

:Default: ``False``






Processors:
 * ``bool_validator`` from `Pydantic`




Frag Ai Iterations
########################################################################






:YAML key: frag_ai_iterations


:Parser key: frag_ai_iterations
:Type: ``int``

:Default: ``False``






Processors:
 * ``int_validator`` from `Pydantic`




Chain Core
########################################################################






:YAML key: chain_core


:Parser key: chain_core
:Type: ``str``

:Default: ``False``






Processors:
 * ``str_validator`` from `Pydantic`




Frag Restart
########################################################################






:YAML key: frag_restart


:Parser key: frag_restart
:Type: ``bool``

:Default: ``False``






Processors:
 * ``bool_validator`` from `Pydantic`




Frag Criteria
########################################################################






:YAML key: frag_criteria


:Parser key: frag_criteria
:Type: ``str``

:Default: ``False``






Processors:
 * ``str_validator`` from `Pydantic`




Frag Output Folder
########################################################################






:YAML key: frag_output_folder


:Parser key: frag_output_folder
:Type: ``str``

:Default: ``False``






Processors:
 * ``str_validator`` from `Pydantic`




Frag Cluster Folder
########################################################################






:YAML key: frag_cluster_folder


:Parser key: frag_cluster_folder
:Type: ``str``

:Default: ``False``






Processors:
 * ``str_validator`` from `Pydantic`




Frag Library
########################################################################






:YAML key: frag_library


:Parser key: frag_library
:Type: ``str``






Processors:
 * ``str_validator`` from `Pydantic`




Frag Core Atom
########################################################################






:YAML key: frag_core_atom


:Parser key: frag_core_atom
:Type: ``str``






Processors:
 * ``str_validator`` from `Pydantic`




Analysis To Point
########################################################################






:YAML key: analysis_to_point


:Parser key: analysis_to_point
:Type: ``float``









*******************************************************
Constraints
*******************************************************


Ca Constr
########################################################################






:YAML key: ca_constr


:Parser key: ca_constr
:Type: ``int``




If value is ``None``,
it will fall back to ``ca_constr`` software setting.
If there is no software setting,
it will default to ``5``.




Processors:
 * ``int_validator`` from `Pydantic`




Ca Interval
########################################################################






:YAML key: ca_interval


:Parser key: ca_interval
:Type: ``float``




If value is ``None``,
it will fall back to ``ca_interval`` software setting.
If there is no software setting,
it will default to ``5``.




Processors:
 * ``float_validator`` from `Pydantic`




Constrain Core Spring
########################################################################






:YAML key: constrain_core_spring


:Parser key: constrain_core_spring
:Type: ``int``

:Default: ``50.0``






Processors:
 * ``int_validator`` from `Pydantic`




External Constraints
########################################################################






:YAML key: external_constraints


:Parser key: external_constraints
:Type: ``str``

:Default: ``[]``








Permissive Metal Constr
########################################################################






:YAML key: permissive_metal_constr


:Parser key: permissive_metal_constr
:Type: ``Any``

:Default: ``list()``








Constrain All Metals
########################################################################






:YAML key: constrain_all_metals


:Parser key: constrain_all_metals
:Type: ``Any``

:Default: ``False``








No Metal Constraints
########################################################################






:YAML key: no_metal_constraints


:Parser key: no_metal_constraints
:Type: ``Any``

:Default: ``False``









*******************************************************
Analysis
*******************************************************


Only Analysis
########################################################################






:YAML key: only_analysis


:Parser key: only_analysis
:Type: ``bool``

:Default: ``False``






Processors:
 * ``bool_validator`` from `Pydantic`





*******************************************************
Metals
*******************************************************


Polarize Metals
########################################################################






:YAML key: polarize_metals


:Parser key: polarize_metals
:Type: ``bool``

:Default: ``False``






Processors:
 * ``bool_validator`` from `Pydantic`




Polarization Factor
########################################################################






:YAML key: polarization_factor


:Parser key: polarization_factor
:Type: ``float``

:Default: ``2.0``






Processors:
 * ``float_validator`` from `Pydantic`





*******************************************************
Custom workflows
*******************************************************


Workflow
########################################################################






:YAML key: workflow


:Parser key: workflow
:Type: ``Any``








Distance
########################################################################






:YAML key: distance


:Parser key: distance
:Type: ``float``






Processors:
 * ``float_validator`` from `Pydantic`





*******************************************************
PPI
*******************************************************


N Components
########################################################################






:YAML key: n_components


:Parser key: n_components
:Type: ``int``


:Tests value: ``3``



If value is  falsy,
it will fall back to ``n_components`` software setting.
If there is no software setting,
it will default to ``10``.




Processors:
 * ``int_validator`` from `Pydantic`




Ppi
########################################################################






:YAML key: ppi


:Parser key: ppi
:Type: ``bool``






Processors:
 * ``bool_validator`` from `Pydantic`




Center Of Interface
########################################################################






:YAML key: center_of_interface


:Parser key: center_of_interface
:Type: ``str``






Processors:
 * ``str_validator`` from `Pydantic`

 * ``validate_atom_string`` from `YamlParserModel`




Protein
########################################################################






:YAML key: protein


:Parser key: protein
:Type: ``str``






Processors:
 * ``str_validator`` from `Pydantic`




Ligand Pdb
########################################################################






:YAML key: ligand_pdb


:Parser key: ligand_pdb
:Type: ``str``






Processors:
 * ``str_validator`` from `Pydantic`




Skip Refinement
########################################################################






:YAML key: skip_refinement


:Parser key: skip_refinement
:Type: ``bool``

:Default: ``False``






Processors:
 * ``bool_validator`` from `Pydantic`





*******************************************************

*******************************************************


N Components
########################################################################






:YAML key: n_components


:Parser key: n_components
:Type: ``int``


:Tests value: ``3``



If value is  falsy,
it will fall back to ``n_components`` software setting.
If there is no software setting,
it will default to ``10``.




Processors:
 * ``int_validator`` from `Pydantic`




Skip Refinement
########################################################################






:YAML key: skip_refinement


:Parser key: skip_refinement
:Type: ``bool``

:Default: ``False``






Processors:
 * ``bool_validator`` from `Pydantic`




Allosteric
########################################################################






:YAML key: allosteric


:Parser key: allosteric
:Type: ``bool``






Processors:
 * ``bool_validator`` from `Pydantic`





*******************************************************
RNA
*******************************************************


Rna
########################################################################






:YAML key: rna


:Parser key: rna
:Type: ``bool``






Processors:
 * ``bool_validator`` from `Pydantic`





*******************************************************
GPCR
*******************************************************


Gpcr Orth
########################################################################






:YAML key: gpcr_orth


:Parser key: gpcr_orth
:Type: ``bool``






Processors:
 * ``bool_validator`` from `Pydantic`




Orthosteric Site
########################################################################






:YAML key: orthosteric_site


:Parser key: orthosteric_site
:Type: ``str``






Processors:
 * ``str_validator`` from `Pydantic`

 * ``validate_atom_string`` from `YamlParserModel`




Initial Site
########################################################################






:YAML key: initial_site


:Parser key: initial_site
:Type: ``str``






Processors:
 * ``str_validator`` from `Pydantic`

 * ``validate_atom_string`` from `YamlParserModel`




Final Site
########################################################################






:YAML key: final_site


:Parser key: final_site
:Type: ``str``






Processors:
 * ``str_validator`` from `Pydantic`

 * ``validate_atom_string`` from `YamlParserModel`




