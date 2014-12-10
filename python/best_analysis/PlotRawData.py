#!/usr/bin/env python
"""Plot a species from the data we get from BEST"""

__author__     = "Nils Wenzler"
__copyright__  = "Nils Wenzler"
__maintainer__ = "Nils Wenzler"
__email__      = "wenzlern@ethz.ch"

### Imports
import pdb
import numpy as np
from PlotHelper import *

def PlotRawData(data, field, species=None):
    # Get a figure
    fig = plt.figure() 
    ax = fig.add_subplot(111)

    PlotFieldOverX(ax, data, args.field, species = args.species)

    fig.colorbar(GetColorBar(ax))

    return fig, ax


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import argparse
    from DataHelper import *

    def CommandLineArgs():
        ### Parse command line arguments
        parser = argparse.ArgumentParser()
        ## Command line arguments
        parser.add_argument('filename', nargs='+', type=str,  
                            help='The mat/npz file with the data')
        
        parser.add_argument('-s', '--species', nargs='*', type=str,
                            help='Provide a species to Filter for: Electrolyte, Cathode...')
        
        parser.add_argument('-f', '--field', nargs='?', type=str,
                            help='Provide a field to plot: [concentration, potential, ...]')
        
        return parser.parse_args()

    args = CommandLineArgs()

    data = ReadData(args.filename[0])

    fig, ax = PlotRawData(data, args.field, species = args.species)
    plt.show()

