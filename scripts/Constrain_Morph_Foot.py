"""
Name: Constrain_Morph_Foot.py

Purpose:
    This script constrains the nodes on the top edge of the morph section of the ankle to the nodes on the bottom edge of
    the shank. This ensures a watertight model.

Usage:
    abaqus cae -nogui Constrain_Morph_Foot.py

Dependencies:
    - numpy - Used for numerical computations

Author:
    Noah Palermo (npalermo@mines.edu)

Last Modified:
    August 23, 2024

Notes:
    Ensure that this script is in the same parent directory as the CAE database.
    
"""

from abaqus import *
from abaqusConstants import *
import numpy as np
import regionToolset
openMdb(
    pathName='foot_model_v16-2021.cae')

#   Initialize abaqus variables
  
a = mdb.models['ankle_model'].rootAssembly
p = mdb.models['ankle_model'].parts['brace']


setNameM = 'Bottom Ring'

morphSet = a.sets[setNameM]
morphElements = morphSet.elements

setNameF = 'Top Ring'
footSet = a.sets[setNameF]
footElements = footSet.elements

fullfootElements = a.instances['smooth_foot_STL'].elements
fullmorphElements = a.instances['smooth_anklemorph_STL'].elements

#   get_centroid returns the centroid of a triangular mesh element given the coordinates of its nodes
def get_centroid(element):
    nodecoords = []
    
    #   Get nodes from element
    nodes = element.getNodes()
    
    #   Iterate through nodes and append coordinates of each to nodecoords
    for i, el in enumerate(nodes):
        coordtuple = nodes[i].coordinates
        nodecoords.append([coordtuple[0], coordtuple[1], coordtuple[2]])
        
    nodecoords = np.array(nodecoords)

    #   Calculate centroid of element
    centroid = np.mean(nodecoords, axis=0)
    
    return centroid
 
#   pdist returns the linear distance between any two points in Cartesian space
def pdist(point1, point2):
    return np.linalg.norm(point1 - point2)
        

morphCentroids = []
footCentroids = []


#   Iterate through every element in the top ring of the foot
for element in footElements:

    #   Get the centroid of the current element
    centroidel = get_centroid(element)
    
    #   Create a surface for the current element
    elementnum = element.label - 30000
    sideElements = fullfootElements[elementnum:elementnum + 1]
    namestr = "Foot" + str(element.label)
    a.Surface(side12Elements=sideElements, name=namestr)
    
    #   Append surface name and centroid of the current element to footCentroids    
    footCentroids.append([namestr, centroidel[0], centroidel[1], centroidel[2]])


#   Iterate through every element in the bottom ring of the morph
for element in morphElements:

    #   Get the centroid of the current element
    centroidel = get_centroid(element)
    
    #   Create a surface for the current element
    elementnum = element.label - 40000
    sideElements = fullmorphElements[elementnum:elementnum + 1]
    namestr = "Morph" + str(element.label)
    a.Surface(side12Elements=sideElements, name=namestr)
    
    #   Append surface name and centroid of the current element to morphCentroids    
    morphCentroids.append([namestr, centroidel[0], centroidel[1], centroidel[2]])
    
    
pairs = []


#   Iterate through every element in morphCentroids
for morphelement in morphCentroids:

        #   Initialize minimum distance found
        mindist = 1e9
        
        #   Get coordinates of the centroid of the current morph element
        morphcoords = np.array([morphelement[1], morphelement[2], morphelement[3]])
        
        #   Iterate through every element in footCentroids
        for footelement in footCentroids:
            
            #   Get coordinates of the centroid of the current foot element
            footcoords = np.array([footelement[1], footelement[2], footelement[3]])
            
            #   Find distance between two current elements
            currdist = pdist(morphcoords, footcoords)
            
            #   Check to see if distance is the closest so far
            if (currdist < mindist):
                mindist = currdist
                closest = footelement[0]
                
        #   Append the two closest elements to a list, along with their linear distance
        pairs.append([closest, morphelement[0], mindist])


session.viewports['Viewport: 1'].assemblyDisplay.setValues(interactions=ON, 
    constraints=ON, connectors=ON, engineeringFeatures=ON, 
    optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)


#   Iterate through every pair in list pairs
for pair in pairs:

    #   Define a tie constraint between each pair
    region1 = a.surfaces[pair[0]]
    region2 = a.surfaces[pair[1]]
    namestr = str(pair[0]) + str(pair[1])
    mdb.models['ankle_model'].Tie(name=namestr, master=region1, slave=region2, 
        positionToleranceMethod=SPECIFIED, positionTolerance=1.0, adjust=OFF, 
        tieRotations=ON, thickness=ON)

#   Save model
mdb.save()