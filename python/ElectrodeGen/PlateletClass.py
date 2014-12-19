#!/usr/bin/env python
"""Stores all relevant info about one platelet"""

__author__     = "Nils Wenzler"
__copyright__  = "Nils Wenzler"
__maintainer__ = "Nils Wenzler"
__email__      = "wenzlern@ethz.ch"

# Imports
import lxml.etree
import lxml.builder
import collections
import random
import numpy as np
import math

def Rotate3DVector(vector, axis, theta):
    vector = np.asarray(vector)
    axis = np.asarray(axis)
    theta = np.asarray(theta)*np.pi/180

    axis = axis/math.sqrt(np.dot(axis, axis))

    a = math.cos(theta/2)
    b, c, d = -axis*math.sin(theta/2)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    RotMat = np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                      [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                      [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

    return np.dot(RotMat, vector)

def GetAxis(direction, angle):
    """Converts a direction vector (3-dim) to 2 axis for geodict"""
    #Direction/Axis1 is perpendicular to the plane of the platelet
    #Perpendiculat/Axis2 is along the first Ray
    axis1 = direction/np.linalg.norm(direction) 
    axis2 = Rotate3DVector([0,-1,0], direction, angle)

    return axis1.tolist(), axis2.tolist()

def GetRayLengths(longdiam, shortdiam, sizedisp=None, shapedisp=None):
    """Creates the raylength and angles from input param"""
    Lengths = [0,0,0,0]

    if shapedisp:
        # Lets keep them symetrical, since we can set the other thing with the angles
        Lengths[0] = random.uniform(longdiam-longdiam*sizedisp,
                                    longdiam+longdiam*sizedisp)
        Lengths[1] = random.uniform(shortdiam-shortdiam*sizedisp,
                                    shortdiam+shortdiam*sizedisp)
        Lengths[2] = random.uniform(longdiam-longdiam*sizedisp,
                                    longdiam+longdiam*sizedisp)
        Lengths[3] = random.uniform(shortdiam-shortdiam*sizedisp,
                                    shortdiam+shortdiam*sizedisp)
        return Lengths

    elif sizedisp:
        # Lets keep them symetrical, since we can set the other thing with the angles
        Lengths[0] = random.uniform(longdiam-longdiam*sizedisp,
                                    longdiam+longdiam*sizedisp)
        Lengths[1] = random.uniform(shortdiam-shortdiam*sizedisp,
                                    shortdiam+shortdiam*sizedisp)
        Lengths[2] = Lengths[0]
        Lengths[3] = Lengths[1]
        return Lengths

    else:
        Lengths[0] = longdiam/2
        Lengths[1] = shortdiam/1
        Lengths[2] = longdiam/2
        Lengths[3] = shortdiam/1
        return Lengths


# Classes
class Platelet:
    """This class saves all the information for one Platelet"""
    def __init__(self, number, position, direction, angle, 
                 longdiam=4e-6, shortdiam=2e-6, thickness=1e-6,
                 sizedisp=None, thicknessdisp=None, shapedisp=None):
        self.Param = collections.OrderedDict()

        #Fill in the info
        self.Name = 'Object' + str(number)
        self.Param['Color'] = random.randint(1,2**4-1)

        self.Param['Type'] = 'PlanarPolyhedron'
        self.Param['Position'] = position
        self.Param['Thickness'] = thickness

        axis1, axis2 = GetAxis(direction, angle)
        self.Param['Axis1'] = axis1
        self.Param['Axis2'] = axis2

        RayLengths = GetRayLengths(longdiam, shortdiam, sizedisp, shapedisp)
        self.Param['Rays'] = 4
        self.Param['RayLength1'] = RayLengths[0]
        self.Param['RayLength2'] = RayLengths[1]
        self.Param['RayLength3'] = RayLengths[2]
        self.Param['RayLength4'] = RayLengths[3]

        self.Param['RayAngle1'] = 90
        self.Param['RayAngle2'] = 90
        self.Param['RayAngle3'] = 90

    def GetPos(self):
        return np.asarray(self.Param['Position'])

    def GetAx1(self):
        return np.asarray(self.Param['Axis1'])

    def GetAx2(self):
        return np.asarray(self.Param['Axis2'])

    
    def SetNr(self, number):
        self.Name = 'Object' + str(number)
    
    def ToString(self):
            # Return a string with dictionary key and value
            ret = '<' + self.Name + '>\n'
            for k, v in self.Param.items():
                ret += str(k) + ' ' + str(v).strip('[] ') + '\n' 
            ret += '</' + self.Name + '>\n'
            return ret
