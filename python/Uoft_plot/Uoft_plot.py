#!/usr/bin/env python
"""Uoft.log file plotter
This script reads the provided Uoft.log file and plots SOC vs Potential"""

__author__     = "Nils Wenzler"
__copyright__  = "Nils Wenzler"
__version__    = '1'
__maintainer__ = "Nils Wenzler"
__email__      = "wenzlern@ethz.ch"

#imports
import argparse
import sys
import matplotlib.pyplot as plt
import numpy as np
from PlotClass import *

## Parse command line arguments
parser = argparse.ArgumentParser()
# Command line arguments
parser.add_argument('filenames', nargs='+', type=str,  
                    help='The folder that contains the Uoft.log file')
parser.add_argument('-l', '--legend', nargs='*', type=str,
                    help='Provide a list of names: [labelA,labelB]')
parser.add_argument('-c', '--crate', nargs='*', type=int,
                    help='C-Rate of the simulation')
parser.add_argument('-t', '--title', nargs='*', type=str,
                    help='Title for the graph')
args = parser.parse_args()

#Set legend to folder names if not specified
legend = args.legend if args.legend else args.filenames

# Define empty data list
time     = []
transf_c = []
cell_pot = []
curr     = []
SOC      = []

# Load the 4 columns starting with the 2nd line
for i in range(len(args.filenames)):
    data = np.loadtxt(args.filenames[i] + '/Uoft.log')
    time.append(np.asarray(data[:,0][:]))
    transf_c.append(np.asarray(data[:,1][:]))
    cell_pot.append(np.asarray(data[:,2][:]))
    curr.append(np.asarray(data[:,3][:]))
    SOC.append((time[i])/(3600/args.crate[i]))


# Define arrays for the individual values

for i in range(len(args.filenames)):
    plt.plot(SOC[i], cell_pot[i], label=legend[i])

plt.grid(True)
#plt.axis([0,1,0,1.5])
plt.xlabel('SOC[-]')
plt.ylabel('Cell Voltage[V]')
plt.legend(loc=2)
plt.title(' '.join(args.title))
plt.show()
