"""
Name: Extract_Rotations.py

Purpose:
    This script extracts the final rotations about the subtalar and talocrural axes from an ODB after analysis. Results
    are printed to the most recent .rpy file.

Usage:
    abaqus cae -nogui Extract_Rotations.py

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

#   Open ODB
o3 = session.openOdb(
    name='sprain.odb')
odb = session.odbs['sprain.odb']

#   Extract talocrural rotation data
xy1 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Connector element relative rotation: CUR1 at Element 70005 in ELSET TALOCRURAL_AXIS', 
    suppressQuery=True, __linkedVpName__='Viewport: 1')
c1 = session.Curve(xyData=xy1)
xyp = session.XYPlot('XYPlot-1')
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
chart.setValues(curvesToPlot=(c1, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)

#   Print talocrural rotation data to .rpy
print('Talocrural Axis rotation:')
print(xy1)
print()
print('Final Rotation:')
print(xy1[-1][-1])
print()

#   Extract subtalar rotation data
xy1 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Connector element relative rotation: CUR1 at Element 70006 in ELSET SUBTALAR_AXIS', 
    suppressQuery=True, __linkedVpName__='Viewport: 1')
c1 = session.Curve(xyData=xy1)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
chart.setValues(curvesToPlot=(c1, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)

#   Print subtalar rotation data to .rpy
print('Subtalar Axis rotation:')
print(xy1)
print()
print('Final Rotation:')
print(xy1[-1][-1])
