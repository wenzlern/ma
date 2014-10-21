#!/usr/bin/env python
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

class _DataItem():
    """This class contains all necessary information to plot a data object"""
    def __init__(self, name, data):
        self.name= name
        self.data = data

class GraphItem():
    """Contains all info for generating one graph"""
    #Define colors and linestyles to iterate through
    colors = ['b','g','r','c','m','y','k','w']
    linestyles = ['', ' ', '--', '-.', '-', ':']

    def __init__(self, ax, xlabel, ylabel, title, grid=True, prange=[], plots=[]):
        self.ax = ax
        self.plots = plots
        self.xlabel= xlabel
        self.ylabel= ylabel
        self.title = title
        self.grid  = grid
        self.prange = prange

    def addplot(self, name, data):
        self.plots.append(_DataItem(name, data))

    def plot(self):
        self.ax.grid(self.grid)
        self.ax.set_xlabel(self.xlabel)
        self.ax.set_ylabel(self.ylabel)
        self.ax.set_title(self.title)

        if self.prange:
            self.ax.axis(self.prange)
            
        for i in range(len(self.plots)):
            for j in range(len(self.plots[i].data)):
                out = self.ax.plot(self.plots[i].data[0][j], 
                                       self.plots[i].data[1][j], 
                                       color=_GraphItem.colors[i], 
                                       linestyle=_GraphItem.linestyles[j], 
                                       label=self.plots[i].name)
        return out

#class Plot():
#    """Create a figure object and fill it in. Provide plot method"""
#
#    def __init__(self):
#        self.fig = plt.figure()
#        self.plot = self.fig.add_subplot(111)
#        self.graphs = []
#
#    def addGraph(self, title, xlabel, ylabel, grid=True, range=[]): 
#        self.graphs.append(_GraphItem(self.plot, xlabel, ylabel, title, grid, range))
#
#    def addData(self, name, data, graphNr=0):
#        self.graphs[graphNr].addplot(name, data)
#
#    def show(self):
#        self.graphs[0].plot()
#        self.fig.show()
#
