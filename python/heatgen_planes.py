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
#data = np.load(args.filename[0])

#The current-density we just have in the form of vectors, 
#get the norm for each cell and scale it by the cell volume
data['AbsCurrent'] = AbsCurrent(data, (0.65e-6)**3)

#Ugly Hack since there is a mistake in the vtkreader that i can only fix tmrw.
#abscurrhack = []
#for i in range(56):
#    abscurrhack.append(data['AbsCurrent'])
#data['AbsCurrent'] = np.asarray(abscurrhack)

#Species to search for
species = ['Cathode', 'Electrolyte']
fields = ['Potential', 'IndividualOCV', 'AbsCurrent']
FilteredData = {}
#Filter for a species, and find heatgen rate
for i in range(len(species)):
    FilteredData.update((FilterForSpecies(data, fields, species[i])))

#Generate some arrays to plot:
MeanNrNeighborsC = []
MeanNrNeighborsE = []

HeatGen = data['AbsCurrent']*(data['IndividualOCV'] 
          - data['Potential'])

#Plot colors
#We need as many colors as xsize
xsize = np.shape(data['Concentration'])[1]
timesteps = np.shape(data['Concentration'])[0]
ColorMap = plt.get_cmap('autumn_r')

plt.subplot(121)
plt.title('Active Material')
plt.grid(True)
plt.ylabel('Mean heat generation')
plt.xlabel('Distance from Current Collector')
plt.gca().set_color_cycle([ColorMap(1.*i/timesteps) for i in range(timesteps)])


plt.subplot(122)
plt.title('Electrolyte')
plt.grid(True)
plt.ylabel('Mean heat generation')
plt.xlabel('Distance from Current Collector')
plt.gca().set_color_cycle([ColorMap(1.*i/timesteps) for i in range(timesteps)])

CathodeIndexMat = FindSpecies(data, 'Cathode')
ElectrolyteIndexMat = FindSpecies(data, 'Electrolyte')

for j in range(timesteps):
    #We have xsize planes 
    TempC = []
    TempE = []    
    for i in range(xsize):
       
        TempC.append(-np.mean(HeatGen[j][i][CathodeIndexMat[i]])) 
        TempE.append(-np.mean(HeatGen[j][i][ElectrolyteIndexMat[i]])) 

    MeanNrNeighborsC.append(np.asarray(TempC))
    MeanNrNeighborsE.append(np.asarray(TempE))
    
    if j > 0:
        plt.subplot(121)
        plt.plot(range(xsize), MeanNrNeighborsC[j], color = ColorMap(1.0*j/timesteps))
        plt.subplot(122)
        plt.plot(range(xsize), MeanNrNeighborsE[j], color = ColorMap(1.0*j/timesteps))

#    if j > 0:
#        if  not(np.allclose(MeanNrNeighborsC[j-1], MeanNrNeighborsC[j], rtol=0.0001)):
#            plt.subplot(211)
#            plt.title('Active Material')
#            plt.plot(range(9)[::-1], MeanNrNeighborsC[j], color = ColorMap(1.*j/xsize))
#        if  not(np.allclose(MeanNrNeighborsE[j-1], MeanNrNeighborsE[j], rtol=0.00001)):
#            plt.subplot(212)
#            plt.title('Electrolyte')
#            plt.plot(range(9)[::-1], MeanNrNeighborsE[j], color = ColorMap(1.*j/xsize))

plt.show()
