"""
Name: Set_Axis_Stiffness.py

Purpose:
    This script is designed define non-linear stiffness curves for the talocrural and subtalar axes of the model. As
    materials are not used for the foot, shank, or ankle, axis-based stiffness is used to approximate passive joint
    stiffness. This is a pre-processing script - included CAE database has been pre-configured, and this script may
    be used to make adjustments to the stiffness curves as desired.

Usage:
    abaqus cae -nogui Set_Axis_Stiffness.py

Dependencies:
    - numpy - Used for numerical computations

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
import numpy as np

#   Open CAE database
executeOnCaeStartup()

openMdb(
    pathName='foot_model_v16-2021.cae')

#   Initialize Abaqus variables
p = mdb.models['ankle_model'].parts['brace']
a = mdb.models['ankle_model'].rootAssembly

#   Set stiffness table resolution. This defines the number of entries in the non-linear stiffness table for both axes
#   used to approximate the stress functions.
resolution = 100

#   Set ROM and define angle vector for talocrural
start = -60
end = 30
phi_a_t = np.linspace(start, end, resolution)

#   Set ROM and define angle vector for subtalar
start = -13
end = 35
phi_a_s = np.linspace(start, end, resolution)

#   Empty lists used to store stiffness data
table_talocrural = []
table_subtalar = []

#   Define start and endpoints as analogues to asymptotic behavior, append start point to table lists
point1t = (10000000, 61*np.pi/180)
point2t = (-10000000, -31*np.pi/180)
point1s = (1000000, (14)*np.pi/180)
point2s = (-1000000, (-26)*np.pi/180)
table_talocrural.append(point1t)
table_subtalar.append(point1s)


#   Iterate through entire ROM for talocrural axis
for element in phi_a_t:

    #   Convert angle to radians
    angledeg = element
    anglerad = element*np.pi/180
    
    #   Calculate moment using stiffness approximation function
    M = -((-0.03)*np.exp(2.8 + 0.125*angledeg) + 0.05*np.exp(-2 - 0.14*angledeg) + 0.47)*1000
    
    #   Create tuple, this is the data type used by Abaqus for stiffness table
    tup = (M, anglerad)
    
    #   Append stiffness and rotation value to table_talocrural
    table_talocrural.append(tup)

#   Iterate through entire ROM for subtalar axis
for element in phi_a_s:

    #   Convert angle to radians
    angledeg = element
    anglerad = element*np.pi/180
    
    #   Calculate moment using stiffness approximation function
    M = -(-0.02*np.exp(-2.5 + 0.35*angledeg) + 0.12*np.exp(-0.7 - 0.38*angledeg) - 0.11)*1000
    
    #   Create tuple, this is the data type used by Abaqus for stiffness table
    tup = (M, anglerad)
    
    #   Append stiffness and rotation value to table_stubtalar
    table_subtalar.append(tup)


#   Append end point to tables
table_talocrural.append(point2t)
table_subtalar.append(point2s)

#   Convert table lists to tuples
table_talocrural = tuple(table_talocrural)
table_subtalar = tuple(table_subtalar)

#   Set table_talocrural as nonlinear stiffness in Abaqus for talocrural connector
elastic_0 = connectorBehavior.ConnectorElasticity(components=(4, ), 
    behavior=NONLINEAR, table=table_talocrural)
mdb.models['ankle_model'].sections['hinge_talocrural'].setValues(
    behaviorOptions = (elastic_0, ) )
mdb.models['ankle_model'].sections['hinge_talocrural'].behaviorOptions[0].ConnectorOptions(
    )

#   Set table_subtalar as nonlinear stiffness in Abaqus for subtalar connector
elastic_1 = connectorBehavior.ConnectorElasticity(components=(4, ), 
    behavior=NONLINEAR, table=table_subtalar)
mdb.models['ankle_model'].sections['hinge_subtalar'].setValues(
    behaviorOptions = (elastic_1, ) )
mdb.models['ankle_model'].sections['hinge_subtalar'].behaviorOptions[0].ConnectorOptions(
    )

#   Save model
mdb.save()
