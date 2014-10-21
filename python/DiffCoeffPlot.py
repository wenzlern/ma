#!/usr/bin/env python
"""Plot the diff coeff from paper doi: 10.1149/1.1872737"""

__author__     = "Nils Wenzler"
__copyright__  = "Nils Wenzler"
__version__    = '1'
__maintainer__ = "Nils Wenzler"
__email__      = "wenzlern@ethz.ch"

#imports
import matplotlib.pyplot as plt
from numpy import *

## Parse command line arguments
#parser = argparse.ArgumentParser()
## Command line arguments
#
#
#parser.parse_args()

# Define arrays for the individual values
Temp = linspace(200,380,num=400)
DiffCoeff = 10**(-4.43-0.22-54/(Temp-(229-5)))

plt.plot(Temp, DiffCoeff)
plt.yscale('log')
plt.grid(True)
plt.xlabel('Temp [K]')
plt.ylabel('D(T) [cm^2/s')
plt.show()
