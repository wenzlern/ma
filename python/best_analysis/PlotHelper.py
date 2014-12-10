"""Plotting functions 
This class provides objects and methods to simplify plotting"""

__author__     = "Nils Wenzler"
__copyright__  = "Nils Wenzler"
__maintainer__ = "Nils Wenzler"
__email__      = "wenzlern@ethz.ch"

### Imports
import pdb
import matplotlib.pyplot as plt
import numpy as np
from DataHelper import *

### HelperFunctions
def ConfigurePlot(ax, xlabel, ylabel, title, legend=None, grid=True):
    """Configures an empty subplot with label, title etc."""
    #Set parameters of the plot
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(grid)

    #Set legend if available
    if legend:
        ax.legend()

def AddPlotData(ax, xdata, ydata, legend=None, color=None):
    """Adds the specified data to an existing plot"""
    if color:
        ax.plot(xdata, ydata, label = legend, color = color)
    else:
        ax.plot(xdata, ydata, label = legend)

def PlotData(ax, xdata, ydata, xlabel, ylabel, title, legend=None, grid=True):
    """Returns a plot object of the data"""
    #Get an empty plot and add data
    ConfigurePlot(ax, xlabel, ylabel, title, legend=legend, grid=grid)
    AddPlotData(ax, xdata, ydata, legend = legend)

def PlotSeries(ax, xdata, ydata, xlabel, ylabel, title, grid=True, legend=None, 
               colormap = 'autumn_r'):
    """Returns a plot object of a series of data"""
    # First get a new plot
    ConfigurePlot(ax, xlabel, ylabel, title, legend=legend, grid=grid)
    #Define colormap
    ColorMap = plt.get_cmap(colormap)
    # We need to know the number of steps
    steps = np.shape(ydata)[0]

    #Set color cycle according to our data
    ax.set_color_cycle([ColorMap(1.*i/steps) for i in range(steps)])

    #Set legend if available or generate one if True
    if legend:
        if type(legend) == type(list()):
            for i in range(steps):
                AddPlotData(ax, xdata, ydata[i], legend=legend[i], 
                            color = ColorMap(1.0*i/steps))
        elif legend == True:
            for i in range(steps):
                AddPlotData(ax, xdata, ydata[i], legend = 't=' + str(i), 
                            color = ColorMap(1.0*i/steps))
    else:
        for i in range(steps):
            AddPlotData(ax, xdata, ydata[i], color = ColorMap(1.0*i/steps))


def GetColorBar(ax, colormap='autumn_r', min=0, max=1):
    cm = plt.cm.ScalarMappable(cmap=plt.get_cmap(colormap), norm=plt.normalize(vmin=min, vmax=max))
    # fake up the array of the scalar mappable...
    cm._A = []
    return cm


def PlotFieldOverX(ax, data, field, species=None, timestep=None):
    #Get the length of x
    xsize = np.shape(data['Concentration'])[1]

    #First we filter for a species if necessary
    if species:
        FiltData = FilterForSpecies(data, field, species)[field][1::]
        #Update field tag
        field = field+species
    else:
        FiltData = data[field][1::]

    #Calculate the mean of each plane, without counting 0s
    #[:-2] repeats the same shape apart from the last 2, which get multiplied
    #Then we can sum over the last index which is now the yz plane
#    PlaneSums = np.sum(FiltData.reshape(FiltData.shape[:-2] +
#                       (np.prod(FiltData.shape[-2:]),)), axis = -1)
#    NrNotZero = np.sum((FiltData!=0).reshape(FiltData.shape[:-2] + 
#                       (np.prod(FiltData.shape[-2:]),)), axis = -1)

    if species:
        title = 'Mean yz plane ' + field + ' in the ' + species + ' as a function of the '
        title += 'distance to current collector of the Anode'
    else:
        title = 'Mean yz plane  ' + field + ' as a function of the ' 
        title += 'distance to current collector of the Anode' 

    #Calculate the mean
#    PlaneMeans = PlaneSums/NrNotZero
    PlaneMeans = CalcYZPlaneMean(FiltData)
    PlotSeries(ax, range(xsize), PlaneMeans,
               'Distance from current collector of Anode', field, title)

 
