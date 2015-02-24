#!/usr/bin/env python
"""This file generates the Open-Circuit-Potential for LFP,
   according to DOI: 10.1149/1.3567007"""

__author__     = "Nils Wenzler"
__copyright__  = "Nils Wenzler"
__maintainer__ = "Nils Wenzler"
__email__      = "wenzlern@ethz.ch"

# Imports
from ParameterFile import *
from math import exp
import sys
import matplotlib.pyplot as plt
import numpy as np
import argparse

def LFPFitting(soc):

    ocp = 3.4323
    ocp += -(0.8428*exp(-80.2493*(1-soc)**1.3198))
    ocp += -(3.2474*1e-6*exp(20.2645*(1-soc)**3.8003))
    ocp += 3.2482*1e-6*exp(20.2646*(1-soc)**3.7995)

    return ocp

def LFPFittingInv(soc):

    ocp = 3.4323
    ocp += -(0.8428*exp(-80.2493*(soc)**1.3198))
    ocp += -(3.2474*1e-6*exp(20.2645*(soc)**3.8003))
    ocp += 3.2482*1e-6*exp(20.2646*(soc)**3.7995)

    return ocp

def LFPFittingInvLower(soc):

    ocp = 3.4323
    ocp += -(0.8428*exp(-80.2493*(soc)**1.3198))
    ocp += -(3.2474*1e-6*exp(20.2645*(soc)**3.8003))
    ocp += 3.2482*1e-6*exp(20.2646*(soc)**3.7995)

    return ocp - 1.5

def LFPFittingTest(soc):

    ocp = 3.45 - 0.2*soc - 3.45*(soc-0.5)**3

    return ocp

def movingaverage(interval, window_size):
    window = np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'same')

## Parse command line arguments
parser = argparse.ArgumentParser()
# Command line arguments
parser.add_argument('filename', nargs='?', type=str,  
                    help='The folder that contains the data')

parser.add_argument('--s', '--start', nargs='?', type=float,  
                    help='The start voltage of the curve')

parser.add_argument('--e', '--end', nargs='?', type=float,  
                    help='The end voltage of the curve')

parser.add_argument('-p', '--plot', action='store_true',
                    help='Keep updating plot')

args = parser.parse_args()


comment = 'LFP SignatureCurve @25 degree celsius, from DOI: 10.1149/1.3567007'

Paramfile = Parametrization('LFP_OCP_RT', comment = comment)

if args.filename:
    data = np.loadtxt(args.filename, skiprows=2)
    var = np.asarray(data[:,1][:])

else:
    Paramfile.AddFunc(LFPFitting, 0, 1, 0.01)
    var = Paramfile.Values


# To calculate the numerical differentiation we need one more value at each
# end of the range.
#var =np.concatenate(([4.2], Paramfile.Values))
#var = movingaverage(var, 4)
step = 1/float(len(var))
#var_c = np.diff(np.hstack(args.voltage,var))/step
var_c = (var[1:] - var[0:-1])/step
print step, var_c[0]


Paramfile.AddValues(var)
Paramfile.AddValues_c(np.hstack((var_c[1],var_c[1:],var_c[-1])))

Paramfile.SaveToFile('LFP_OCP_RT.c')
    
if args.plot:
    plt.plot(range(len(Paramfile.Values)), Paramfile.Values)
    plt.plot(range(len(Paramfile.Values_c)), Paramfile.Values_c)
    plt.show()
