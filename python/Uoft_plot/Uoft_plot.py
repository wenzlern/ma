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
from numpy import *
from PlotClass import *

## Parse command line arguments
#parser = argparse.ArgumentParser()
## Command line arguments
#
#
#parser.parse_args()


# Load the 4 columns starting with the 2nd line
data = loadtxt(sys.argv[1], skiprows = 1)

# Define arrays for the individual values
time, transf_c, cell_pot, curr = data[:,0], data[:,1], data[:,2], data[0,3]
SOC = (transf_c)/curr

Plot = Plot()
Plot.addGraph('SOC vs Cell Voltage', 'SOC[-]', 'Cell Voltage [V]')
Plot.addData('298K', [SOC,cell_pot])
Plot.show()

#plt.plot(SOC, cell_pot)
#plt.grid(True)
#plt.axis([0,1,0,1.5])
#plt.xlabel('SOC[-]')
#plt.ylabel('Cell Voltage[V]')
#plt.show()
