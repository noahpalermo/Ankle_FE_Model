"""
Name: Set_Axis_Sensors.py

Purpose:
    This script adds sensors to the talocrural and subtalar axes in Abaqus.

Usage:
    abaqus cae -nogui Set_Axis_Sensors.py

Dependencies:
    - none

Author:
    Noah Palermo (npalermo@mines.edu)

Last Modified:
    August 25, 2024

Notes:
    Ensure that this script is in the same parent directory as the CAE database.
    
"""

from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup

#   Open CAE database
executeOnCaeStartup()

    openMdb('foot_model_v16-2021.cae')

#   Initialize Abaqus variables
p = mdb.models['ankle_model'].parts['brace']
a = mdb.models['ankle_model'].rootAssembly

#   Insert subtalar sensor
regionDef=mdb.models['ankle_model'].rootAssembly.sets['subtalar_axis']
mdb.models['ankle_model'].HistoryOutputRequest(name='Subtalar_Rotation', 
    createStepName='defined_pose', variables=('UR', ), region=regionDef, 
    sectionPoints=DEFAULT, rebar=EXCLUDE, sensor=ON)

#   Insert talocrural sensor
regionDef=mdb.models['ankle_model'].rootAssembly.sets['talocrural_axis']
mdb.models['ankle_model'].HistoryOutputRequest(name='Talocrural_Rotation', 
    createStepName='defined_pose', variables=('UR', ), region=regionDef, 
    sectionPoints=DEFAULT, rebar=EXCLUDE, sensor=ON)

#   Save model
mdb.save()