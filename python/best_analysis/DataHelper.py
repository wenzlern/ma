#!/usr/bin/env python
"""Some helper functions to analyse the best output"""

__author__     = "Nils Wenzler"
__copyright__  = "Nils Wenzler"
__maintainer__ = "Nils Wenzler"
__email__      = "wenzlern@ethz.ch"

#imports
import numpy as np
#import scipy.io as spio

Materials = { 'Anode' : 0,
              'Cathode' : 100,
              'Electrolyte' : 200,
              'CCAnode' : 300,
              'CCCathode' : 400,
}

def ReadData(file):
    if file.split('.')[-1] == 'mat':
        return spio.loadmat(file)
    else:
        data = np.load(file)
        # Due to the way numpy stores we need to sanitize dicts
        res = dict()
        res.update(data)
        res['Parameters'] = res['Parameters'][()]
        res['ParamApplication'] = res['ParamApplication'][()]
        res['InterfaceAreas'] = res['InterfaceAreas'][()]

        return res


def FilterForSpecies(data, fields, species):
    """Returns a dict of of the specified fields with all 
       but the specified species set to 0. The dict key of the output is
       'Field'."""
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
        output[fields[j]+species] = res

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
    for x in range(1,xr-1):
        for y in range(1, yr-1):
            for z in range(1, zr-1):
                mat = data['MaterialIdentifier'][0, x, y, z]
               # cube = data['MaterialIdentifier'][0, x-1:x+2, y-1:y+2, z-1:z+2]
               # res[x,y,z] = np.sum(cube == cube[1,1,1])
                res[x,y,z] += 1.0*(mat == data['MaterialIdentifier'][0, x+1, y, z])
                res[x,y,z] += 1.0*(mat == data['MaterialIdentifier'][0, x-1, y, z])
                res[x,y,z] += 1.0*(mat == data['MaterialIdentifier'][0, x, y+1, z])
                res[x,y,z] += 1.0*(mat == data['MaterialIdentifier'][0, x, y-1, z])
                res[x,y,z] += 1.0*(mat == data['MaterialIdentifier'][0, x, y, z+1])
                res[x,y,z] += 1.0*(mat == data['MaterialIdentifier'][0, x, y, z-1])
    return res

def AbsCurrent(data, area = 1):
    """Returns a matrix of scalars of the absolute of the vectorized
       current density. If volume is specified the Current density is
       converted to an absolute current"""
    #First reduce dimension to 2-dim matrix for easier math
    cellarray = np.reshape(data['CurrentDensity'],
                          (np.prod(np.shape(data['CurrentDensity'])[0:-1]),3))
    #The new version of numpy.linalg.norm has an axis property and is faster
    #NormArray = (np.norm(cellarray, axis=1))*volume
    #We have an old numpy/python version, oh well:
    NormArray = (np.sum(np.abs(cellarray)**2,axis=1)**0.5)*area

    #Convert back to original form minus the last dimension (reduced in norm)
    res = np.reshape(NormArray, np.shape(data['CurrentDensity'])[0:-1])

    return res


def CalcElectrodeMean(data):
    """Calculate the mean of the whole electrode, without counting 0s"""
    #[:-3] repeats the same shape apart from the last 3, which get multiplied
    #Then we can sum over the last index which is now the yz plane
    PlaneSums = np.sum(data.reshape(data.shape[:-3] +
                       (np.prod(data.shape[-3:]),)), axis = -1)
    NrNotZero = np.sum((data!=0).reshape(data.shape[:-3] + 
                       (np.prod(data.shape[-3:]),)), axis = -1)

    return PlaneSums/NrNotZero
     

def CalcYZPlaneMean(data):
    """Calculate the mean of each plane, without counting 0s"""
    #[:-2] repeats the same shape apart from the last 2, which get multiplied
    #Then we can sum over the last index which is now the yz plane
    PlaneSums = np.sum(data.reshape(data.shape[:-2] +
                       (np.prod(data.shape[-2:]),)), axis = -1)
    NrNotZero = np.sum((data!=0).reshape(data.shape[:-2] + 
                       (np.prod(data.shape[-2:]),)), axis = -1)

    return PlaneSums/NrNotZero
    
 
    
