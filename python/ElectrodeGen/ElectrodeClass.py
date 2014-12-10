#!/usr/bin/env python
"""Stores all relevant info about one Electrode"""

__author__     = "Nils Wenzler"
__copyright__  = "Nils Wenzler"
__maintainer__ = "Nils Wenzler"
__email__      = "wenzlern@ethz.ch"

# Imports
import lxml.etree
import lxml.builder
import collections
import random
from PlateletClass import Platelet

# Classes
class Electrode:
    """This class saves all the information for a full Electrode description"""
    def __init__(self, size, voxelsize=5e-07):
        self.Param = collections.OrderedDict()
        self.Objects = []

        #Fill in the World info
        self.Param['PeriodicX'] = False
        self.Param['PeriodicY'] = False
        self.Param['PeriodicZ'] = False

        self.Param['OriginX'] = 0
        self.Param['OriginY'] = 0
        self.Param['OriginZ'] = 0

        self.Param['VoxelLength'] = voxelsize
        self.Param['WorldMode'] = 'VoxelNumber'

        self.Param['NX'] =  size[0]
        self.Param['NY'] =  size[1]
        self.Param['NZ'] =  size[2]

        self.Param['OverlapMode'] = 'OverlapColor'
        self.Param['Color'] = 0

    def ToString(self):
        # Header
        ret = '<GadAddDialog>\nUnit 2\nKeepStructure 0\n'
        
        # Add the World info
        ret += '<World>\n'
        # Return a string with dictionary key and value
        for k, v in self.Param.items():
            ret += str(k) + ' ' + str(v) + '\n' 
        ret += '</World>\n'

        #Objects count
        ret += 'NumberOfObjects ' + str(len(self.Objects)) + '\n'

        # Add the Objects
        for s in self.Objects:
            ret += s.ToString()

        # Add end Header
        ret += '</GadAddDialog>'
        return ret

    def SaveToFile(self,Filename):
        f = open(Filename, 'w')
        f.write(self.ToString())
        f.close()
