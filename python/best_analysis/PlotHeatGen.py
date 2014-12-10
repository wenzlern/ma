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
from CalculateHeatGen import *

def NormalizeVector(v):
    norm = np.linalg.norm(v)
    return v/norm

def PlotCellVoltage(ax, data):
    ax.set_ylabel('V')
    AddPlotData(ax, range(len(data['Time'])), data['CellPot'], 'CellVoltage')

def PlotOverPotential(ax, data, species=None):
    if species:
        temp = (FilterForSpecies(data, ['Overpotential'], species))['Overpotential'+species]
    else:
        temp = data['Overpotential']

    MeanOverPotential = CalcElectrodeMean(temp)

    AddPlotData(ax, range(len(MeanOverPotential)), MeanOverPotential, 'Overpotential'+species)

def PlotPotential(ax, data, species=None):
    if species:
        temp = (FilterForSpecies(data, ['Potential'], species))['Potential'+species]
    else:
        temp = data['Potential']

    MeanPotential = CalcElectrodeMean(temp)

    AddPlotData(ax, range(len(MeanPotential)), MeanPotential, 'Potential'+species)

def PlotCurrent(ax, data, species=None):
    data['AbsCurrent'] = AbsCurrent(data, area = (0.65e-4)**2)
    if species:
        temp = (FilterForSpecies(data, ['AbsCurrent'], species))['AbsCurrent'+species]
    else:
        temp = data['AbsCurrent']

    MeanCurrent = CalcElectrodeMean(temp)

    AddPlotData(ax, range(len(MeanCurrent)), MeanCurrent, 'Current'+species)
    return MeanCurrent

def PlotJouleHeating(data, species=None):
    # Get a figure
    fig = plt.figure() 
    ax = fig.add_subplot(111)

    JouleHeating = SumJouleHeating(data)
    
    PlotData(ax, range(len(JouleHeating)), JouleHeating, 'Timesteps', 'W/cm2',
             'Joule Heating over Time')

    return fig, ax

def PlotElectroChemHeating(data):
    # Get a figure
    fig = plt.figure() 
    ax = fig.add_subplot(111)

    ElectroChemHeating = SumElectroChemHeating(data, 'Cathode')
    
    PlotData(ax, range(len(ElectroChemHeating)), ElectroChemHeating, 'Timesteps', 'W/cm2',
             'Electro-Chemical Heating over Time')

    return fig, ax
    
def PlotAllHeating(data):
    fig = plt.figure() 
    ax = fig.add_subplot(111)

    JouleHeatingElectrolyte = SumJouleHeating(data, species='Electrolyte')
    JouleHeatingCathode = SumJouleHeating(data, species='Cathode')
    JouleHeating = JouleHeatingCathode + JouleHeatingElectrolyte
    ElectroChemHeating = SumElectroChemHeating(data, 'Cathode', area=float(data['InterfaceAreas']['Cathode']))

    # Combine the data to put it into a series plot
    ydata = np.asarray([JouleHeatingElectrolyte, JouleHeatingCathode, JouleHeating,
                        ElectroChemHeating, JouleHeating+ElectroChemHeating])
    print np.shape(ydata)
    
    PlotSeries(ax, range(len(JouleHeating)), 
               ydata, 'Timesteps [30s]', 'W/cm2', 'Heating over Time', legend=['JouleElectrolyte',
               'JouleCathode', 'JouleTotal', 'EChem', 'Sum'], colormap = 'jet')

    return fig, ax


def PlotAllAbsHeating(data):
    fig = plt.figure() 
    ax = fig.add_subplot(111)

    JouleHeatingElectrolyte = SumJouleHeating(data, area=(0.65e-4)**2, species='Electrolyte')
    JouleHeatingCathode = SumJouleHeating(data, area=(0.65e-4)**2, species='Cathode')
    JouleHeating = JouleHeatingCathode + JouleHeatingElectrolyte
    ElectroChemHeating = SumElectroChemHeating(data, 'Cathode', area=float(data['InterfaceAreas']['Cathode']))

    # Combine the data to put it into a series plot
    ydata = np.asarray([JouleHeatingElectrolyte, JouleHeatingCathode, JouleHeating,
                        ElectroChemHeating, JouleHeating+ElectroChemHeating])
    print np.shape(ydata)
    
    PlotSeries(ax, range(len(JouleHeating)), 
               ydata, 'Timesteps [30s]', 'W/cm2', 'Heating over Time', legend=['JouleElectrolyte',
               'JouleCathode', 'JouleTotal', 'EChem', 'Sum'], colormap = 'jet')

    return fig, ax

def PlotDifferentialSig(data):
    fig = plt.figure() 
    ax = fig.add_subplot(111)

    JouleHeatingElectrolyte = SumJouleHeating(data, area=(0.65e-4)**2, species='Electrolyte')
    JouleHeatingCathode = SumJouleHeating(data, area=(0.65e-4)**2, species='Cathode')
    ElectroChemHeating = SumElectroChemHeating(data, 'Cathode', area=float(data['InterfaceAreas']['Cathode']))
    DiffSig = JouleHeatingCathode - JouleHeatingElectrolyte

    # Combine the data to put it into a series plot
    ydata = np.asarray([JouleHeatingElectrolyte, JouleHeatingCathode, DiffSig, ElectroChemHeating])

    PlotSeries(ax, range(len(DiffSig)), 
               ydata, 'Timesteps [30s]', 'W/cm2', 'Heat generation rate over Time', legend=['JouleElectrolyte',
               'JouleCathode', 'DiffSig', 'EChem'], colormap = 'jet')

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
        
        parser.add_argument('-t', '--type', nargs='*', type=str,
                            help='List of the types of heating to plot,' + 
                            'default is All: {All, Joule, ElectroChem}')
        
        return parser.parse_args()

    args = CommandLineArgs()

    data = ReadData(args.filename[0])

    #fig, ax = PlotJouleHeating(data, species = args.species)
    #fig, ax = PlotElectroChemHeating(data)
    #fig, ax = PlotAllHeating(data)
    fig, ax = PlotAllAbsHeating(data)
    #fig, ax = PlotDifferentialSig(data)


    ax2 = ax.twinx()
    ax2.set_ylabel('V')
    PlotCurrent(ax, data, 'Cathode')
    PlotCurrent(ax, data, 'Electrolyte')

    PlotPotential(ax2, data, 'Cathode')
    PlotPotential(ax2, data, 'Electrolyte')
    
    ax.set_yscale('log')

    #Get all plotted lines and labels
    lines, labels = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()

    ax2.legend(lines + lines2, labels + labels2, loc=0)


    plt.show()

