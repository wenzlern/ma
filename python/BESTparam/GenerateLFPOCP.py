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

comment = 'LFP SignatureCurve @25 degree celsius, from DOI: 10.1149/1.3567007'

Paramfile = Parametrization('LFP_OCP_RT', comment = comment)

Paramfile.AddFunc(LFPFittingInv, 0, 1, 0.01)
# To calculate the numerical differentiation we need one more value at each
# end of the range.
var =np.concatenate(([2.5], Paramfile.Values, [LFPFittingInv(1.01)]))
var_c = var[2:len(var)] - var[0:-2]
#Paramfile.AddValues_c(np.diff(np.concatenate(([2.5], Paramfile.Values, [LFPFittingInv(1.01)])))[1:-1])
Paramfile.AddValues_c(var_c)

Paramfile.SaveToFile('LFP_OCP_RT.c')

if sys.argv[1] == '-p':
    print Paramfile.Values[0], Paramfile.Values[1]
    print Paramfile.Values_c[0], Paramfile.Values_c[1]
    plt.plot(range(len(Paramfile.Values)), Paramfile.Values)
    plt.plot(range(len(Paramfile.Values_c)), Paramfile.Values_c)
    plt.show()
