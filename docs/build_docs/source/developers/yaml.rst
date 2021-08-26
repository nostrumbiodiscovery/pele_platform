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

:Default: ``''``






Processors:
 * ``str_validator`` from `Pydantic`

 * ``validate_residue_name`` from `YamlParserModel`

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




Restart
########################################################################






:YAML key: restart


:Parser key: restart
:Type: ``bool``






Processors:
 * ``bool_validator`` from `Pydantic`




Singularity Exec
########################################################################






:YAML key: singularity_exec


:Parser key: singularity_exec
:Type: ``str``






Processors:
 * ``str_validator`` from `Pydantic`




Inter Step Logger
########################################################################






:YAML key: inter_step_logger


:Parser key: inter_step_logger
:Type: ``bool``






Processors:
 * ``bool_validator`` from `Pydantic`

 * ``set_interstep_logger`` from `SimulationParamsModel`





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




Templates
########################################################################






:YAML key: templates


:Parser key: templates
:Type: ``Any``

:Default: ``'/home/agruzka/pele_platform/pele_platform/PeleTemplates'``








Solvent Template
########################################################################






:YAML key: solvent_template


:Parser key: solvent_template
:Type: ``str``






Processors:
 * ``str_validator`` from `Pydantic`




Ext Rotamers
########################################################################






:YAML key: None, value from `rotamers` parser


:Parser key: ext_rotamers
:Type: ``str``








Mpi Params
########################################################################






:YAML key: mpi_params


:Parser key: mpi_params
:Type: ``str``






Processors:
 * ``str_validator`` from `Pydantic`

 * ``format_mpi_params`` from `SimulationParamsModel`




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




Exit Clust / Clust
########################################################################




.. warning::
    This is a candidate for deprecation.



:YAML key: exit_clust


:Parser key: clust
:Type: ``int``


:Tests value: ``2``





Processors:
 * ``int_validator`` from `Pydantic`




One Exit
########################################################################




.. warning::
    This is a candidate for deprecation.



:YAML key: one_exit


:Parser key: one_exit
:Type: ``Any``








Box Metric
########################################################################




.. warning::
    This is a candidate for deprecation.



:YAML key: box_metric


:Parser key: box_metric
:Type: ``Any``








Time
########################################################################




.. warning::
    This is a candidate for deprecation.



:YAML key: time


:Parser key: time
:Type: ``Any``








Nosasa
########################################################################




.. warning::
    This is a candidate for deprecation.



:YAML key: nosasa


:Parser key: nosasa
:Type: ``Any``








Perc Sasa
########################################################################




.. warning::
    This is a candidate for deprecation.



:YAML key: perc_sasa


:Parser key: perc_sasa
:Type: ``Any``








Seed
########################################################################






:YAML key: seed


:Parser key: seed
:Type: ``int``

:Default: ``generate_random_seed()``


:Tests value: ``12345``





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




.. warning::
    This is a candidate for deprecation.



:YAML key: log


:Parser key: log
:Type: ``Any``








Nonrenum
########################################################################




.. warning::
    This is a candidate for deprecation.



:YAML key: nonrenum


:Parser key: nonrenum
:Type: ``Any``








Pele Exec
########################################################################






:YAML key: pele_exec


:Parser key: pele_exec
:Type: ``str``

:Default: ``'/scratch/PELE/V1.7.1/bin/Pele_mpi'``






Processors:
 * ``str_validator`` from `Pydantic`




Pele Data
########################################################################






:YAML key: pele_data


:Parser key: pele_data
:Type: ``str``

:Default: ``'/scratch/PELE/V1.7.1/Data'``






Processors:
 * ``str_validator`` from `Pydantic`




Pele Documents
########################################################################






:YAML key: pele_documents


:Parser key: pele_documents
:Type: ``str``

:Default: ``'/scratch/PELE/V1.7.1/Documents'``






Processors:
 * ``str_validator`` from `Pydantic`




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
:Type: ``int``




If value is  falsy,
it will fall back to ``anm_modes_change`` software setting.
If there is no software setting,
it will default to ``4``.




Processors:
 * ``int_validator`` from `Pydantic`




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
:Type: ``Union[List[str], str]``








Analyse
########################################################################






:YAML key: analyse


:Parser key: analyse
:Type: ``bool``

:Default: ``True``






Processors:
 * ``bool_validator`` from `Pydantic`




Mae
########################################################################






:YAML key: mae


:Parser key: mae
:Type: ``bool``

:Default: ``False``






Processors:
 * ``bool_validator`` from `Pydantic`




Spawning Condition
########################################################################

For min or maximising epsilon




:YAML key: spawning_condition


:Parser key: spawning_condition
:Type: ``str``




If value is  falsy,
it will fall back to ``spawning_condition`` software setting.




Processors:
 * ``str_validator`` from `Pydantic`

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




Pele License
########################################################################






:YAML key: pele_license


:Parser key: pele_license
:Type: ``str``

:Default: ``'/scratch/PELE/V1.7.1/licenses'``






Processors:
 * ``str_validator`` from `Pydantic`




License
########################################################################






:YAML key: None, value from `pele_license` parser


:Parser key: license
:Type: ``str``






Processors:
 * ``str_validator`` from `Pydantic`




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




Interaction Restrictions
########################################################################






:YAML key: interaction_restrictions


:Parser key: interaction_restrictions
:Type: ``dict``






Processors:
 * ``parse_interaction_restrictions`` from `SimulationParamsModel`




Terminal Constr
########################################################################






:YAML key: terminal_constr


:Parser key: terminal_constr
:Type: ``float``




If value is  falsy,
it will fall back to ``terminal_constr`` software setting.
If there is no software setting,
it will default to ``5.0``.




Processors:
 * ``float_validator`` from `Pydantic`

 * ``assert_positive_integer`` from `YamlParserModel`




Minimum Steps
########################################################################






:YAML key: minimum_steps


:Parser key: minimum_steps
:Type: ``bool``






Processors:
 * ``bool_validator`` from `Pydantic`

 * ``set_minimum_steps`` from `SimulationParamsModel`




Site Finder Global
########################################################################






:YAML key: site_finder_global


:Parser key: site_finder_global
:Type: ``bool``






Processors:
 * ``bool_validator`` from `Pydantic`




Site Finder Local
########################################################################






:YAML key: site_finder_local


:Parser key: site_finder_local
:Type: ``bool``






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
:Type: ``str``






Processors:
 * ``str_validator`` from `Pydantic`

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




External Templates
########################################################################





:YAML key: None, value calculated in simulation params

:Parser key: external_templates
:Type: ``Any``




If value is  falsy,
it will fall back to ``templates`` software setting.
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

:Default: ``'/opt/schrodinger2020-1/utilities/python'``






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








Sidechain Perturbation
########################################################################





:YAML key: None, value calculated in simulation params

:Parser key: sidechain_perturbation
:Type: ``str``






Processors:
 * ``str_validator`` from `Pydantic`

 * ``only_with_sidechain_perturbation`` from `SimulationParamsModel`




Met Interaction Restrictions
########################################################################





:YAML key: None, value calculated in simulation params

:Parser key: met_interaction_restrictions
:Type: ``str``

:Default: ``''``






Processors:
 * ``str_validator`` from `Pydantic`




Covalent Sasa
########################################################################





:YAML key: None, value calculated in simulation params

:Parser key: covalent_sasa
:Type: ``str``






Processors:
 * ``str_validator`` from `Pydantic`

 * ``only_with_sidechain_perturbation`` from `SimulationParamsModel`




Max Trials For One
########################################################################





:YAML key: None, value calculated in simulation params

:Parser key: max_trials_for_one
:Type: ``Any``








Conformation Perturbation
########################################################################





:YAML key: None, value calculated in simulation params

:Parser key: conformation_perturbation
:Type: ``str``






Processors:
 * ``str_validator`` from `Pydantic`

 * ``set_conformation_perturbation`` from `SimulationParamsModel`




Equilibration Mode
########################################################################





:YAML key: None, value calculated in simulation params

:Parser key: equilibration_mode
:Type: ``str``






Processors:
 * ``str_validator`` from `Pydantic`




Water Ids To Track
########################################################################





:YAML key: None, value calculated in simulation params

:Parser key: water_ids_to_track
:Type: ``str``

:Default: ``[]``








Inputs Dir
########################################################################





:YAML key: None, value calculated in simulation params

:Parser key: inputs_dir
:Type: ``str``






Processors:
 * ``str_validator`` from `Pydantic`




Residue Type
########################################################################





:YAML key: None, value calculated in simulation params

:Parser key: residue_type
:Type: ``str``






Processors:
 * ``str_validator`` from `Pydantic`





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

 * ``check_divisibility`` from `YamlParserModel`




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
:Type: ``Union[List[float], str]``




If value is  falsy,
it will fall back to ``cluster_values`` software setting.
If there is no software setting,
it will default to ``'[1.75, 2.5, 4, 6]'``.






Cluster Conditions
########################################################################






:YAML key: cluster_conditions


:Parser key: cluster_conditions
:Type: ``Union[List[float], str]``




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

 * ``check_divisibility`` from `YamlParserModel`




Core
########################################################################






:YAML key: core


:Parser key: core
:Type: ``Union[int, List[str]]``






Processors:
 * ``set_core_constraints`` from `SimulationParamsModel`




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

Maximum number of flexible side chains in the ligand.




:YAML key: n


:Parser key: n
:Type: ``int``

:Default: ``10000``






Processors:
 * ``int_validator`` from `Pydantic`




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






Charge Parametrization Method
########################################################################






:YAML key: charge_parametrization_method


:Parser key: charge_parametrization_method
:Type: ``str``






Processors:
 * ``str_validator`` from `Pydantic`




Exclude Terminal Rotamers
########################################################################






:YAML key: exclude_terminal_rotamers


:Parser key: exclude_terminal_rotamers
:Type: ``bool``

:Default: ``True``






Processors:
 * ``bool_validator`` from `Pydantic`




Use Peleffy
########################################################################






:YAML key: use_peleffy


:Parser key: use_peleffy
:Type: ``bool``

:Default: ``False``






Processors:
 * ``bool_validator`` from `Pydantic`





*******************************************************
Protein preparation
*******************************************************


Skip Preprocess / No Ppp
########################################################################




.. warning::
    This is a candidate for deprecation.



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




Nonstandard
########################################################################






:YAML key: nonstandard


:Parser key: nonstandard
:Type: ``str``




If value is  falsy,
it will fall back to ``nonstandard`` software setting.
If there is no software setting,
it will default to ``[]``.







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








Box Type
########################################################################






:YAML key: box_type


:Parser key: box_type
:Type: ``str``






Processors:
 * ``str_validator`` from `Pydantic`





*******************************************************
Metrics
*******************************************************


Rmsd Pdb / Native
########################################################################






:YAML key: rmsd_pdb


:Parser key: native
:Type: ``str``

:Default: ``''``






Processors:
 * ``str_validator`` from `Pydantic`




Atom Dist
########################################################################






:YAML key: atom_dist


:Parser key: atom_dist
:Type: ``str``

:Default: ``list()``






Processors:
 * ``validate_atom_string`` from `YamlParserModel`





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
MSM
*******************************************************


Msm
########################################################################




.. warning::
    This is a candidate for deprecation.



:YAML key: msm


:Parser key: msm
:Type: ``bool``






Processors:
 * ``bool_validator`` from `Pydantic`




Lagtime
########################################################################




.. warning::
    This is a candidate for deprecation.



:YAML key: lagtime


:Parser key: lagtime
:Type: ``Any``








Msm Clust
########################################################################




.. warning::
    This is a candidate for deprecation.



:YAML key: msm_clust


:Parser key: msm_clust
:Type: ``Any``









*******************************************************
Rescoring
*******************************************************


Rescoring
########################################################################






:YAML key: rescoring


:Parser key: rescoring
:Type: ``bool``






Processors:
 * ``bool_validator`` from `Pydantic`





*******************************************************
In Out
*******************************************************


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
:Type: ``bool``






Processors:
 * ``bool_validator`` from `Pydantic`




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
:Type: ``float``




If value is  falsy,
it will fall back to ``water_temp`` software setting.
If there is no software setting,
it will default to ``5000``.




Processors:
 * ``float_validator`` from `Pydantic`




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
:Type: ``float``




If value is  falsy,
it will fall back to ``water_constr`` software setting.
If there is no software setting,
it will default to ``0``.




Processors:
 * ``float_validator`` from `Pydantic`




Water Trials
########################################################################






:YAML key: water_trials


:Parser key: water_trials
:Type: ``int``




If value is  falsy,
it will fall back to ``water_trials`` software setting.
If there is no software setting,
it will default to ``10000``.




Processors:
 * ``int_validator`` from `Pydantic`




Water Radius
########################################################################






:YAML key: water_radius


:Parser key: water_radius
:Type: ``float``

:Default: ``6.0``






Processors:
 * ``float_validator`` from `Pydantic`




Water Empty Selector
########################################################################






:YAML key: water_empty_selector


:Parser key: water_empty_selector
:Type: ``bool``

:Default: ``False``






Processors:
 * ``bool_validator`` from `Pydantic`




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




.. warning::
    This is a candidate for deprecation.



:YAML key: frag


:Parser key: frag
:Type: ``bool``






Processors:
 * ``bool_validator`` from `Pydantic`




COMligandConstraint / Com
########################################################################






:YAML key: COMligandConstraint


:Parser key: com
:Type: ``float``




If value is  falsy,
it will fall back to ``COMligandConstraint`` software setting.
If there is no software setting,
it will default to ``0``.




Processors:
 * ``float_validator`` from `Pydantic`




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






Processors:
 * ``str_validator`` from `Pydantic`




Frag Ligands
########################################################################






:YAML key: frag_ligands


:Parser key: frag_ligands
:Type: ``Any``








Growing Steps
########################################################################






:YAML key: growing_steps


:Parser key: growing_steps
:Type: ``int``

:Default: ``6``






Processors:
 * ``int_validator`` from `Pydantic`




Steps In Gs / Frag Steps
########################################################################






:YAML key: steps_in_gs


:Parser key: frag_steps
:Type: ``int``

:Default: ``3``






Processors:
 * ``int_validator`` from `Pydantic`




Sampling Steps / Frag Eq Steps
########################################################################






:YAML key: sampling_steps


:Parser key: frag_eq_steps
:Type: ``int``

:Default: ``20``






Processors:
 * ``int_validator`` from `Pydantic`




Protocol
########################################################################






:YAML key: protocol


:Parser key: protocol
:Type: ``str``

:Default: ``''``






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

:Default: ``'L'``






Processors:
 * ``str_validator`` from `Pydantic`




Frag Restart
########################################################################






:YAML key: frag_restart


:Parser key: frag_restart
:Type: ``str``

:Default: ``''``






Processors:
 * ``str_validator`` from `Pydantic`




Frag Criteria
########################################################################






:YAML key: frag_criteria


:Parser key: frag_criteria
:Type: ``str``

:Default: ``'Binding Energy'``






Processors:
 * ``str_validator`` from `Pydantic`




Frag Output Folder
########################################################################






:YAML key: frag_output_folder


:Parser key: frag_output_folder
:Type: ``str``

:Default: ``'growing_steps'``






Processors:
 * ``str_validator`` from `Pydantic`




Frag Cluster Folder
########################################################################






:YAML key: frag_cluster_folder


:Parser key: frag_cluster_folder
:Type: ``str``

:Default: ``'clustering_PDBs'``






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

 * ``validate_frag_core_atom`` from `YamlParserModel`




Analysis To Point
########################################################################






:YAML key: analysis_to_point


:Parser key: analysis_to_point
:Type: ``float``








Fragment Atom
########################################################################






:YAML key: fragment_atom


:Parser key: fragment_atom
:Type: ``str``






Processors:
 * ``str_validator`` from `Pydantic`




Frag Restart Libraries
########################################################################






:YAML key: frag_restart_libraries


:Parser key: frag_restart_libraries
:Type: ``bool``






Processors:
 * ``bool_validator`` from `Pydantic`





*******************************************************
Constraints
*******************************************************


Ca Constr
########################################################################






:YAML key: ca_constr


:Parser key: ca_constr
:Type: ``float``




If value is  falsy,
it will fall back to ``ca_constr`` software setting.
If there is no software setting,
it will default to ``0.5``.




Processors:
 * ``float_validator`` from `Pydantic`

 * ``assert_positive_integer`` from `YamlParserModel`




Ca Interval
########################################################################






:YAML key: ca_interval


:Parser key: ca_interval
:Type: ``int``




If value is  falsy,
it will fall back to ``ca_interval`` software setting.
If there is no software setting,
it will default to ``10``.




Processors:
 * ``int_validator`` from `Pydantic`




Constrain Core
########################################################################

String of SMILES or SMARTS to constrain.




:YAML key: constrain_core


:Parser key: constrain_core
:Type: ``str``




If value is  falsy,
it will fall back to ``constrain_core`` software setting.




Processors:
 * ``str_validator`` from `Pydantic`




Constrain Core Spring
########################################################################






:YAML key: constrain_core_spring


:Parser key: constrain_core_spring
:Type: ``float``

:Default: ``50.0``






Processors:
 * ``float_validator`` from `Pydantic`




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
:Type: ``bool``






Processors:
 * ``bool_validator`` from `Pydantic`




Constrain All Metals
########################################################################






:YAML key: constrain_all_metals


:Parser key: constrain_all_metals
:Type: ``bool``

:Default: ``False``






Processors:
 * ``bool_validator`` from `Pydantic`




No Metal Constraints
########################################################################






:YAML key: no_metal_constraints


:Parser key: no_metal_constraints
:Type: ``bool``

:Default: ``False``






Processors:
 * ``bool_validator`` from `Pydantic`




Constraint Level
########################################################################






:YAML key: constraint_level


:Parser key: constraint_level
:Type: ``int``






Processors:
 * ``int_validator`` from `Pydantic`

 * ``parse_constraint_level`` from `YamlParserModel`





*******************************************************
PCA
*******************************************************


Pca
########################################################################

Path to PCA file




:YAML key: pca


:Parser key: pca
:Type: ``str``






Processors:
 * ``str_validator`` from `Pydantic`

 * ``format_pca`` from `SimulationParamsModel`





*******************************************************
Advanced
*******************************************************


Perturbation
########################################################################






:YAML key: perturbation


:Parser key: perturbation
:Type: ``bool``






Processors:
 * ``bool_validator`` from `Pydantic`

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




Top Clusters Criterion
########################################################################






:YAML key: top_clusters_criterion


:Parser key: top_clusters_criterion
:Type: ``str``

:Default: ``'interaction_25_percentile'``






Processors:
 * ``str_validator`` from `Pydantic`

 * ``check_selection_criterion`` from `YamlParserModel`




Kde
########################################################################






:YAML key: kde


:Parser key: kde
:Type: ``bool``






Processors:
 * ``bool_validator`` from `Pydantic`




Kde Structs
########################################################################






:YAML key: kde_structs


:Parser key: kde_structs
:Type: ``int``

:Default: ``1000``






Processors:
 * ``int_validator`` from `Pydantic`




Plot Filtering Threshold
########################################################################






:YAML key: plot_filtering_threshold


:Parser key: plot_filtering_threshold
:Type: ``float``






Processors:
 * ``float_validator`` from `Pydantic`




Clustering Filtering Threshold
########################################################################






:YAML key: clustering_filtering_threshold


:Parser key: clustering_filtering_threshold
:Type: ``float``

:Default: ``0.25``






Processors:
 * ``float_validator`` from `Pydantic`




Clustering Method
########################################################################






:YAML key: clustering_method


:Parser key: clustering_method
:Type: ``str``




If value is  falsy,
it will fall back to ``clustering_method`` software setting.
If there is no software setting,
it will default to ``'meanshift'``.




Processors:
 * ``str_validator`` from `Pydantic`




Cluster Representatives Criterion
########################################################################






:YAML key: cluster_representatives_criterion


:Parser key: cluster_representatives_criterion
:Type: ``str``




If value is  falsy,
it will fall back to ``cluster_representatives_criterion`` software setting.
If there is no software setting,
it will default to ``'interaction_5_percentile'``.




Processors:
 * ``str_validator`` from `Pydantic`




Bandwidth
########################################################################






:YAML key: bandwidth


:Parser key: bandwidth
:Type: ``float``




If value is  falsy,
it will fall back to ``bandwidth`` software setting.
If there is no software setting,
it will default to ``2.5``.




Processors:
 * ``float_validator`` from `Pydantic`




Max Top Clusters
########################################################################






:YAML key: max_top_clusters


:Parser key: max_top_clusters
:Type: ``int``

:Default: ``8``






Processors:
 * ``int_validator`` from `Pydantic`




Min Population
########################################################################






:YAML key: min_population


:Parser key: min_population
:Type: ``float``




If value is  falsy,
it will fall back to ``min_population`` software setting.
If there is no software setting,
it will default to ``0.01``.




Processors:
 * ``float_validator`` from `Pydantic`




Max Top Poses
########################################################################






:YAML key: max_top_poses


:Parser key: max_top_poses
:Type: ``int``

:Default: ``100``






Processors:
 * ``int_validator`` from `Pydantic`





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




.. warning::
    This is a candidate for deprecation.



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
Site finder
*******************************************************


N Components
########################################################################




.. warning::
    This is a candidate for deprecation.



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




Site Finder
########################################################################






:YAML key: site_finder


:Parser key: site_finder
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





*******************************************************
Covalent docking
*******************************************************


Covalent Residue
########################################################################






:YAML key: covalent_residue


:Parser key: covalent_residue
:Type: ``str``






Processors:
 * ``str_validator`` from `Pydantic`

 * ``validated_residue_string`` from `YamlParserModel`




Nonbonding Radius
########################################################################






:YAML key: nonbonding_radius


:Parser key: nonbonding_radius
:Type: ``float``






Processors:
 * ``float_validator`` from `Pydantic`




Perturbation Trials
########################################################################






:YAML key: perturbation_trials


:Parser key: perturbation_trials
:Type: ``int``






Processors:
 * ``int_validator`` from `Pydantic`




Refinement Angle
########################################################################






:YAML key: refinement_angle


:Parser key: refinement_angle
:Type: ``float``






Processors:
 * ``float_validator`` from `Pydantic`




Covalent Docking Refinement
########################################################################






:YAML key: covalent_docking_refinement


:Parser key: covalent_docking_refinement
:Type: ``bool``






Processors:
 * ``bool_validator`` from `Pydantic`





*******************************************************
Ligand conformations
*******************************************************


Ligand Conformations
########################################################################






:YAML key: ligand_conformations


:Parser key: ligand_conformations
:Type: ``str``






Processors:
 * ``str_validator`` from `Pydantic`




Conformation Freq
########################################################################






:YAML key: conformation_freq


:Parser key: conformation_freq
:Type: ``int``

:Default: ``4``




If value is  falsy,
it will fall back to ``conformation_freq`` software setting.




Processors:
 * ``int_validator`` from `Pydantic`

 * ``only_with_conformation_perturbation`` from `SimulationParamsModel`




Overlap Factor Conformation
########################################################################






:YAML key: overlap_factor_conformation


:Parser key: overlap_factor_conformation
:Type: ``float``






Processors:
 * ``float_validator`` from `Pydantic`





*******************************************************
Saturated mutagenesis
*******************************************************


Saturated Mutagenesis
########################################################################






:YAML key: saturated_mutagenesis


:Parser key: saturated_mutagenesis
:Type: ``bool``






Processors:
 * ``bool_validator`` from `Pydantic`




Cpus Per Mutation
########################################################################






:YAML key: cpus_per_mutation


:Parser key: cpus_per_mutation
:Type: ``int``


:Tests value: ``2``





Processors:
 * ``int_validator`` from `Pydantic`




