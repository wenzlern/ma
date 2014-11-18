#!/usr/bin/env python
"""Script to calculate heat generation rate in the battery"""

__author__     = "Nils Wenzler"
__copyright__  = "Nils Wenzler"
__maintainer__ = "Nils Wenzler"
__email__      = "wenzlern@ethz.ch"

#imports
import numpy as np
import scipy.io as spio
import matplotlib.pyplot as plt
import argparse
from HelperFunctions import *

### Parse command line arguments
parser = argparse.ArgumentParser()
## Command line arguments
parser.add_argument('filename', nargs='+', type=str,  
                    help='The mat file with the data')

#parser.add_argument('-l', '--legend', nargs='*', type=str,
#                    help='Provide a list of names: labelA labelB ...')
#
#parser.add_argument('-a', '--axis', nargs='*', type=float,
#                    help='Provide axis range: Xmin Xmax Ymin Ymax')
#
#parser.add_argument('--nolegend', action='store_true',
#                    help='Disable the legend')
#
#parser.add_argument('-u', '--update', action='store_true',
#                    help='Keep updating plot')
#
#parser.add_argument('-c', '--crate', nargs='*', type=int,
#                    help='C-Rate of the simulation')
#
#parser.add_argument('-t', '--title', nargs='*', type=str,
#                    help='Title for the graph')
args = parser.parse_args()


#Load file, just temporary a mat...
print args.filename
data = spio.loadmat(args.filename[0])

#Find neighbors and append them to the dictionary
data['NrSameNeighbors'] = FindNrNeighbors(data)

#The current-density we just have in the form of vectors, 
#get the norm for each cell and scale it by the cell volume
data['AbsCurrent'] = AbsCurrent(data, (0.65e-6)**3)

#Ugly Hack since there is a mistake in the vtkreader that i can only fix tmrw.
abscurrhack = []
for i in range(56):
    abscurrhack.append(data['AbsCurrent'])
data['AbsCurrent'] = np.asarray(abscurrhack)

#Species to search for
species = ['Cathode', 'Electrolyte']
fields = ['Potential', 'IndividualOCV', 'AbsCurrent', 'NrSameNeighbors']
FilteredData = {}
#Filter for a species, and find heatgen rate
for i in range(len(species)):
    FilteredData.update((FilterForSpecies(data, fields, species[i])))

#Generate some arrays to plot:
MeanNrNeighborsC = []
MeanNrNeighborsE = []

HeatGen = data['AbsCurrent']*(data['IndividualOCV'] 
          - data['Potential'])
for j in range(56):
    #We have a maximum of 8 neighbors
    TempC = []
    TempE = []    
    for i in range(9):
        CathodeIndexMat = ((FilteredData['NrSameNeighborsCathode']-i) == 0)
        ElectrolyteIndexMat = ((FilteredData['NrSameNeighborsElectrolyte']-i) == 0)
        
        TempC.append(np.mean(HeatGen[j][CathodeIndexMat])) 
        TempE.append(np.mean(HeatGen[j][ElectrolyteIndexMat])) 

    MeanNrNeighborsC.append(np.asarray(TempC))
    MeanNrNeighborsE.append(np.asarray(TempE))

    plt.plot(range(9)[::-1], MeanNrNeighborsC[j], color = str(j/100))
    plt.plot(range(9)[::-1], MeanNrNeighborsE[j], color = str(j/100),linestyle = '--')

plt.grid(True)

plt.ylabel('Mean heat generation')
plt.xlabel('Number of neighbors with same species')

plt.show()
