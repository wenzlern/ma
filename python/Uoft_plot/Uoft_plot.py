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
import scipy.integrate
from scipy.spatial import distance
from PlotClass import *

def spacedmarks(x, y, Nmarks, data_ratio=None):
    if data_ratio is None:
        data_ratio = plt.gca().get_data_ratio()
    
    dydx = gradient(y, x[1])
    dxdx = gradient(x, x[1])*data_ratio
    arclength = scipy.integrate.cumtrapz(sqrt(dydx**2 + dxdx**2), x)
    #We dont have the initial=0 option due to version. Append 0 at beginning
    arclength = np.hstack(([0],arclength))
    marks = linspace(0, max(arclength), Nmarks)
    markx = interp(marks, arclength, x)
    marky = interp(markx, x, y)
    return markx, marky

def FindIdealSpacing(x, y):
    dist = []
    for i in range(len(x)-1):
       dist.append(distance.euclidean([x[i],y[i]],[x[i+1],y[i+1]]))
    return min(dist), max(dist)

## Parse command line arguments
parser = argparse.ArgumentParser()
# Command line arguments
parser.add_argument('filenames', nargs='+', type=str,  
                    help='The folder that contains the Uoft.log file')

parser.add_argument('-l', '--legend', nargs='*', type=str,
                    help='Provide a list of names: labelA labelB ...')

parser.add_argument('-a', '--axis', nargs='*', type=float,
                    help='Provide axis range: Xmin Xmax Ymin Ymax')

parser.add_argument('--nolegend', action='store_true',
                    help='Disable the legend')

parser.add_argument('-u', '--update', action='store_true',
                    help='Keep updating plot')

parser.add_argument('-c', '--crate', nargs='*', type=float,
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
#for i in range(len(args.filenames)):
#    x, y = spacedmarks(SOC[i], cell_pot[i], 30,1.5)
#    plt.plot(x, y, label=legend[i], marker='o')
#

plt.grid(True)
if args.axis:
    plt.axis(args.axis)

plt.xlabel('x in Li_xFePO4[-]')
plt.ylabel('Cell Voltage[V]')

if not args.nolegend:
    plt.legend(loc=2)

plt.title(' '.join(args.title))
plt.show()
