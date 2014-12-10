#!/usr/bin/env python
"""Helper functions to calculate heat generation in batteries"""

__author__     = "Nils Wenzler"
__copyright__  = "Nils Wenzler"
__maintainer__ = "Nils Wenzler"
__email__      = "wenzlern@ethz.ch"

### Imports
import numpy as np
from DataHelper import *

### Calculate
## Joule Heating
def JouleHeating(data, area = 1, species=None):
    """Calculates the Joule heating term of the input data,
       either per area or absolute"""
    # We just filter the Potential, since we multiply [0*x=0]
    if species:
        temp = (FilterForSpecies(data, ['Potential'], species))['Potential'+species]
    else:
        temp = data['Potential']

    return (temp*AbsCurrent(data)*area)

def SumJouleHeating(data, area = 1, species=None):
    """Return the sum of all Joule heating in the specified data"""
    temp = JouleHeating(data, area, species)
    #Flatten to 2-dim and sum over the last dimension
    return np.sum(temp.reshape(temp.shape[0],np.prod(temp.shape[1::])),axis=-1)

## Electrochemical heating
## Here we can be rather lazy at the moment
## If we want to take it from the data we need to find the 
## Voxels with interfaces and the current normal to that face.
## However, the sim-log gives us the total area and we have a fixed
## Exchange current density, so we use this at the moment
def SumElectroChemHeating(data, species, area=1):
    ExchangeCurrDens = float(data['Parameters']['Systems']['MaterialParameters']['ActiveParticle'][species+'1']['BVrate_k'])

    #Get Overpotential of cathode
    data.update(FilterForSpecies(data, ['Overpotential'], species))

    #To get the mean overpotential of all voxels with an interface we need the
    #number of neighbors with same material
    NrNeighbors = FindNrNeighbors(data)
    #Interface are all with less than 8 Neighbors
    InterfaceVoxels = (NrNeighbors < 8)*NrNeighbors

    #ElectrodeInterfaceVox = np.tile(np.shape(data['Overpotential'+species])[0],
    #                        InterfaceVoxels)* data['Overpotential'+species]

    ElectrodeInterfaceVox = InterfaceVoxels*data['Overpotential'+species]

    MeanOverPot = CalcElectrodeMean(ElectrodeInterfaceVox)
    return MeanOverPot*ExchangeCurrDens*area

