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

parser.add_argument('-s', '--species', nargs='*', type=str,
                    help='Provide a species to plot: Electrolyte, Cathode...')

parser.add_argument('-f', '--field', nargs='*', type=str,
                    help='Provide a field to plot: [concentration, potential, ...]')

args = parser.parse_args()


#Load file, just temporary a mat...
print args.filename
data = spio.loadmat(args.filename[0])
#data = np.load(args.filename[0])

#Plot colors and matrix sizes
xsize = np.shape(data['Concentration'])[1]
timesteps = np.shape(data['Concentration'])[0]
ColorMap = plt.get_cmap('autumn_r')


#Species to search for
FilteredData = {}
#Filter for a species
for i in range(len(args.species)):
    FilteredData.update((FilterForSpecies(data, args.field, args.species[i])))

plt.grid(True)
plt.ylabel(args.field)
plt.xlabel('Distance from Current Collector')
plt.gca().set_color_cycle([ColorMap(1.*i/timesteps) for i in range(timesteps)])

CathodeIndexMat = FindSpecies(data, 'Cathode')
ElectrolyteIndexMat = FindSpecies(data, 'Electrolyte')

FieldData = data[args.field[0]]
PlaneMean = []

for j in range(timesteps):
    #We have xsize planes 
    Temp = []
    for i in range(xsize):
        Temp.append(-np.mean(FieldData[j][i][CathodeIndexMat[i]])) 

    PlaneMean.append(np.asarray(Temp))
    plt.plot(range(xsize), PlaneMean[j], color = ColorMap(1.0*j/timesteps))

plt.show()
