#!/usr/bin/env python
"""Some helper functions to analyse the best output"""

__author__     = "Nils Wenzler"
__copyright__  = "Nils Wenzler"
__maintainer__ = "Nils Wenzler"
__email__      = "wenzlern@ethz.ch"

#imports
import numpy as np
#from enum import Enum

Materials = { 'Anode' : 0,
              'Cathode' : 100,
              'Electrolyte' : 200,
              'CCAnode' : 300,
              'CCCathode' : 400,
}

def FilterForSpecies(data, fields, species):
    """Returns a dict of of the specified fields with all 
       but the specified species set to 0. The dict key of the output is
       'FieldSpecies'."""
    output = {}
    for j in range(len(fields)):
        res = np.zeros(np.shape(data[fields[j]]))
        if type(species) == type(str()):
            MaterialID = Materials[species]
        else:
            MaterialID = species

        if len(np.shape(data[fields[j]])) < 4:
            res += (((data['MaterialIdentifier'][0]-MaterialID)==0)*data[fields[j]])
        else:
            res += (((data['MaterialIdentifier']-MaterialID)==0)*data[fields[j]])
        output[fields[j] + species] = res

    return output


def FindSpecies(data, species):
    """Returns a Boolean matrix with True for the fields of that species""" 
    if type(species) == type(str()):
        MaterialID = Materials[species]
    else:
        MaterialID = species

    return ((data['MaterialIdentifier']-MaterialID)==0)[0]


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

def AbsCurrent(data, volume = 1):
    """Returns a matrix of scalars of the absolute of the vectorized
       current density. If volume is specified the Current density is
       converted to an absolute current"""
    #First reduce dimension to 2-dim matrix for easier math
    cellarray = np.reshape(data['CurrentDensity'],
                          (np.prod(np.shape(data['CurrentDensity'])[0:-1]),3))
    #For the absolute value we need to factor in that volume is expected in SI,
    #but out data is (most probably...) in A/cm^3. Therefor volume*10**6
    #The new version of numpy.linalg.norm has an axis property and is faster
    #NormArray = (np.norm(cellarray, axis=1))*volume*10**6
    #We have an old numpy/python version, oh well:
    NormArray = (np.sum(np.abs(cellarray)**2,axis=1)**0.5)*volume*10**6

    #Convert back to original form minus the last dimension (reduced in norm)
    res = np.reshape(NormArray, np.shape(data['CurrentDensity'])[0:-1])

    return res


#def ExtractSpecies(data, field, species):
#    """Returns a numpy array cropped to the size of the specified species."""
    

    
      


