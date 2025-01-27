# Force Controlled Abaqus FE Model

A finite element model of the human ankle aimed at allowing for efficient, iterative design of external orthosis devices.

## Getting Started

The files hosted in this repository are organized in separate directories for clarity. All files in subdirectories must be
placed in the same parend directory for proper usage (e.g. all Python scripts under /Ankle_FE_Model/scripts must be moved to /Ankle_FE_Model.)

### Prerequisites
- Python
- Intel Fortran Compiler
- Abaqus Standard 2021

### Steps
1. Configure Abaqus job
    The included .cae database is pre-configured. The scripts found in /Ankle_FE_Model/scripts may be used to modify the model parameters,
    such as the passive stiffness of anatomical axes, or for postproccesing. Re-write the input file within CAE after any changes have been made.

    Usage of Python Scripts

    The included Python scripts are used to automate processes that can be performed manually in the CAE GUI. They may be submitted at the command line
    in the appropriate directory:

    abaqus cae -nogui Script_Name.py

    More details on usage may be found in the header of each included script.

2. Submitting Abaqus job
    The job must be submitted through Intel Fortran Compiler using the following syntax:

    abaqus -j jobname user=DynamicMorph.f90

## Files
- /cae_files/foot_model_v16-2021.cae: The CAE database of the ankle model.
- /cae_files/foot_model_v16_2021.jnl: The journal file of the CAE database.
- /cae_files/sprain.inp: The input file for the Abaqus job. This will be overwritten should you decide to configure your own job.

- /include_files/morph_include_20dfx_30sup_bcs.inc: Defines the torques applied about the anatomical axes of the model and links the displacement of nodes within the ankle morph to an amplitude computed by DynamicMorph.f90.
- /include_files/MorphData.inc: Defines arrays containing the nodes in the ankle, used by DynamicMorph.f90.

- /scripts/Constrain_Morph_Foot.py: Script for preproccesing - adding tie constraints between the uppermost nodes of the ankle morph and the lowest nodes of the shank.
- /scripts/Extract_Rotations.py: Script for postprocessing - prints the final rotations about the anatomical axes to the most recent .rpy file.
- /scripts/Set_Axis_Sensors.py: Script for preprocessing - defines sensors for the anatomical axes.
- /scripts/Set_Axis_Stiffnes.py: Script for preprocessing - defines stiffness curves for the anatomical axes.

- /subroutines/DynamicMorph.f90: User subroutine allowing for deformation of the ankle morph section of the model.