"""Plotting functions 
This class provides objects and methods to simplify plotting"""

__author__     = "Nils Wenzler"
__copyright__  = "Nils Wenzler"
__version__    = '1'
__maintainer__ = "Nils Wenzler"
__email__      = "wenzlern@ethz.ch"

#imports
import matplotlib.pyplot as plt
from numpy import *

class Plot():
    """Contains all info for generating one graph"""

    #Define colors and linestyles to iterate through
    colors = ['b','g','r','c','m','y','k','w']
    linestyles = ['', ' ', '--', '-.', '-', ':']

    def __init__(self, ax, xlabel, ylabel, titel, xdata, ydata, grid=True)
        

