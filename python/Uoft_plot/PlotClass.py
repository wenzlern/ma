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

class _GraphItem():
    """Contains all info for generating one graph"""
    #Define colors and linestyles to iterate through
    colors = ['b','g','r','c','m','y','k','w']
    linestyles = ['', ' ', '--', '-.', '-', ':']

    def __init__(self, figure, xlabel, ylabel, title, grid=True, range=[], plots=[]):
        self.figure = figure
        self.plots = plots
        self.xlabel= xlabel
        self.ylabel= ylabel
        self.title = title
        self.grid  = grid
        self.range = range

    def addplot(self, name, data):
        self.plots.append(_DataItem(name, data))

    def plot(self):
        self.figure.grid(grid)
        self.figure.xlabel(xlabel)
        self.figure.ylabel(ylabel)
        self.figure.title(title)

        if range:
            self.figure.axis(range)
            
        for i in range(len(plots)):
            for j in range(self.plots[i].data):
                self.figure.plot(self.plots[i].data[1,j], 
                                 self.plots[i].data[2,j], 
                                 color=colors[i], 
                                 linestyle=linestyles[j], 
                                 label=self.plots[i].name)

class Plot():
    """Create a figure object and fill it in. Provide plot method"""

    def __init__(self):
        self.fig = plt.figure()
        self.plot = self.fig.add_subplot(111)
        self.graphs = []

    def addGraph(self, title, xlabel, ylabel, grid=True, range=[]): 
        self.graphs.append(_GraphItem(self.fig, xlabel, ylabel, title, grid, range))

    def addData(self, name, data, graphNr=0):
        self.graphs[graphNr].addplot(name, data)

    def show(self):
        self.graphs[0].plot()
        self.fig.show()

