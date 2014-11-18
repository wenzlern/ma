#!/usr/bin/env python
"""Some helper functions to analyse the best output"""

__author__     = "Nils Wenzler"
__copyright__  = "Nils Wenzler"
__maintainer__ = "Nils Wenzler"
__email__      = "wenzlern@ethz.ch"

#imports
import numpy as np
#from enum import Enum

Materials = { Anode : 0,
              Cathode : 100,
              Electrolyte : 200,
              CCAnode : 300,
              CCCathode : 400,
}

def FilterForSpecies(data, fields, species):
    """Returns a numpy array of the size of the specified field with all 
       but the specified species set to 0. 
       Species can be a list to include several species."""
    output = {}
    for j in range(len(fields)):
        res = np.zeros(np.shape(data[fields[j]]))
        for i in range(len(species)):
            if type(species[i]) is str:
                species[i] = Materials(species[i])

            if len(np.shape(data[fields[j]])) < 4:
                res += (((data['MaterialIdentifier'][0]-species[i])==0)*data[fields[j]])
            else:
                res += (((data['MaterialIdentifier']-species[i])==0)*data[fields[j]])
        output[fields[j]] = res


def FindNrNeighbors(data):
    """Returns a matrix with the same shape of data with the number of neighbors
       with the same material for each voxel."""
    xr, yr, zr = np.shape(data['MaterialIdentifier'][0])

    res = np.zeros((xr, yr, zr))
    cube = np.zeros((3,3,3))
    for x in range(1,xr):
        for y in range(1, yr):
            for z in range(1, zr):
                cube = data['MaterialIdentifier'][0, x-1:x+1, y-1:y+1, z-1:z+1]
                res[x,y,z] = np.sum(cube == cube[1,1,1])

    return res



#def ExtractSpecies(data, field, species):
#    """Returns a numpy array cropped to the size of the specified species."""
    

    
      


