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

comment = 'LFP SignatureCurve @25 degree celsius, from DOI: 10.1149/1.3567007'

Paramfile = Parametrization('LFP_OCP_RT', comment = comment)

Paramfile.AddFunc(LFPFitting, 0, 1, 0.01)
Paramfile.AddValues_c(np.diff(Paramfile.Values))

Paramfile.SaveToFile('LFP_OCP_RT.c')

if sys.argv[1] == '-p':
    plt.plot(range(len(Paramfile.Values)), Paramfile.Values)
    plt.plot(range(len(Paramfile.Values_c)), Paramfile.Values_c)
    plt.show()
