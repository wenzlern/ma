#!/usr/bin/env python
"""Script to calculate heat generation rate in the battery"""

__author__     = "Nils Wenzler"
__copyright__  = "Nils Wenzler"
__maintainer__ = "Nils Wenzler"
__email__      = "wenzlern@ethz.ch"

#imports
import numpy as np
import scipy.io as spio
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
data = spio.loadmat(args.filenames)

#Find neighbors and append them to the dictionary
data['NrSameNeighbors'] = FindNrNeighbors(data)

#Species to search for
species = ['Cathode', 'Electrolyte']
fields = ['']
FilteredData = {}
#Filter for a species, and find heatgen rate
for i in range(len(species)):
    FilteredData.update((FilterForSpecies(data, ))
