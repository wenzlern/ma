#!/usr/bin/env python
"""Stores all relevant info about one platelet"""

__author__     = "Nils Wenzler"
__copyright__  = "Nils Wenzler"
__maintainer__ = "Nils Wenzler"
__email__      = "wenzlern@ethz.ch"

# Imports
#import lxml.etree
#import lxml.builder
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

def normvec(vector):
    """ Returns the unit vector of the vector."""
    return vector / np.linalg.norm(vector)

def angle(v1, v2):

    v1_u = normvec(v1)
    v2_u = normvec(v2)

    angle = np.arccos(np.dot(v1_u, v2_u))
    if np.isnan(angle):
        if (v1_u == v2_u).all():
            return 0.0
        else:
            return np.pi
    return angle

def GetAxis(direction, angle):
    """Converts a direction vector (3-dim) to 2 axis for geodict"""
    #Direction/Axis1 is perpendicular to the plane of the platelet
    #Perpendiculat/Axis2 is along the first Ray
    axis1 = direction/np.linalg.norm(direction) 
    axis2 = Rotate3DVector([0,-1,0], direction, angle)

    return (normvec(axis1)).tolist(), (normvec(axis2)).tolist()

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
        Lengths[0] = random.uniform(longdiam/2-longdiam/2*sizedisp,
                                    longdiam/2+longdiam/2*sizedisp)
        Lengths[1] = random.uniform(shortdiam/2-shortdiam/2*sizedisp,
                                    shortdiam/2+shortdiam/2*sizedisp)
        Lengths[2] = Lengths[0]
        Lengths[3] = Lengths[1]
        return Lengths

    else:
        Lengths[0] = longdiam/2
        Lengths[1] = shortdiam/1
        Lengths[2] = longdiam/2
        Lengths[3] = shortdiam/1
        return Lengths


# Calculate plane from 3 points
def PlaneFrom3Points(p1, p2, p3):
    # Assemble points to matrix
    P = np.transpose(np.vstack((p1, p2, p3)))

    # Calculate the coefficients
    ab = p2-p1
    ac = p3-p1
    A = normvec(np.cross(ab,ac))
    d = -np.dot(A, p1)

    return A, d

def AxbEntryFrom3Points(p1,p2,p3):
    # If the vector points in the negative direction of [1,1,1]
    # we need to multiply A and b with (-1) to flip the sign of
    # that respective Ax<=b entry
    A, d = PlaneFrom3Points(p1, p2, p3)
    if angle(A, np.asarray([1,1,1])) > (np.pi/2):
        A = -A
        d = -d

    return A, d


# Classes
class Platelet:
    """This class saves all the information for one Platelet"""
    def __init__(self, number, position, direction, angle, 
                 longdiam=4e-6, shortdiam=2e-6, thickness=1e-6,
                 sizedisp=None, thicknessdisp=None, shapedisp=None, color=None):
        self.Param = collections.OrderedDict()

        #Fill in the info
        self.Name = 'Object' + str(number)
        if color:
            self.Param['Color'] = random.randint(1,2**4-1)
        else:
            self.Param['Color'] = 1 

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

        self.Modified = False
        
        # Matrix description of the platelet
        self.A = np.zeros((6,3)) 
        self.b = np.zeros((6,1))


    # Calculate the definition of the paltelet via hyperplanes
    def GetPlateletDef(self):
        # We need the 3rd axis, so for easier use we calculate it here
        Axis3 = Rotate3DVector(self.Param['Axis2'], self.Param['Axis1'], 90)
        
        # Define the corners of the platelet clockwise (up and down)
        P11 = np.asarray(self.Param['Position'])
        P11 += (self.Param['RayLength1'] * np.asarray(self.Param['Axis2']))
        P11 += (self.Param['Thickness']  * np.asarray(self.Param['Axis1']))
                                                                     
        P12 = np.asarray(self.Param['Position'])
        P12 += (self.Param['RayLength1'] * np.asarray(self.Param['Axis2']))
        P12 -= (self.Param['Thickness'] *  np.asarray(self.Param['Axis1']))
                                                                     
        P21 = np.asarray(self.Param['Position'])
        P21 += (self.Param['RayLength2'] * Axis3)
        P21 += (self.Param['Thickness'] *  np.asarray(self.Param['Axis1']))
                                                                     
        P22 = np.asarray(self.Param['Position'])
        P22 += (self.Param['RayLength2'] * Axis3)              
        P22 -= (self.Param['Thickness'] *  np.asarray(self.Param['Axis1']))
                                                                     
        P31 = np.asarray(self.Param['Position'])
        P31 -= (self.Param['RayLength3'] * np.asarray(self.Param['Axis2']))
        P31 += (self.Param['Thickness'] *  np.asarray(self.Param['Axis1']))
                                                                     
        P32 = np.asarray(self.Param['Position'])
        P32 -= (self.Param['RayLength3'] * np.asarray(self.Param['Axis2']))
        P32 -= (self.Param['Thickness'] *  np.asarray(self.Param['Axis1']))
                                                                     
        P41 = np.asarray(self.Param['Position'])
        P41 -= (self.Param['RayLength4'] * Axis3)              
        P41 += (self.Param['Thickness'] *  np.asarray(self.Param['Axis1']))
                                                                     
        P42 = np.asarray(self.Param['Position'])
        P42 -= (self.Param['RayLength4'] * Axis3)              
        P42 -= (self.Param['Thickness'] *  np.asarray(self.Param['Axis1']))

        # Now we calculate the hyperplane coefficients
        self.A[0], self.b[0] = AxbEntryFrom3Points(P12, P11, P21)#(P11, P12, P21)
        self.A[1], self.b[1] = AxbEntryFrom3Points(P22, P21, P31)#(P21, P22, P31)
        self.A[2], self.b[2] = AxbEntryFrom3Points(P32, P31, P41)#(P31, P32, P41)
        self.A[3], self.b[3] = AxbEntryFrom3Points(P42, P41, P11)#(P41, P42, P11)
        # Point order changed for A[4] so that the normal vector 
        # points to the center of the platelet
        self.A[4], self.b[4] = AxbEntryFrom3Points(P11, P21, P31)
        self.A[5], self.b[5] = AxbEntryFrom3Points(P12, P22, P32)

#        self.A[0], self.b[0] = PlaneFrom3Points(P12, P11, P21)
#        self.A[1], self.b[1] = PlaneFrom3Points(P22, P21, P31)
#        self.A[2], self.b[2] = PlaneFrom3Points(P32, P31, P41)
#        self.A[3], self.b[3] = PlaneFrom3Points(P42, P41, P11)
#        # Point order changed for A[4] so that the normal vector 
#        # points to the center of the platelet
#        self.A[4], self.b[4] = PlaneFrom3Points(P11, P21, P31)
#        self.A[5], self.b[5] = PlaneFrom3Points(P12, P22, P32)


    def GetPos(self):
        return np.asarray(self.Param['Position'])

    def GetAx1(self):
        return np.asarray(self.Param['Axis1'])

    def GetAx2(self):
        return np.asarray(self.Param['Axis2'])
    
    def SetAx1(self, ax):
        self.Param['Axis1'] = ax
        self.Modified = True

    def SetAx2(self, ax):
        self.Param['Axis2'] = ax
        self.Modified = True

    def SetNr(self, number):
        self.Name = 'Object' + str(number)
        self.Modified = True
    
    def ToString(self):
            # Return a string with dictionary key and value
            ret = '<' + self.Name + '>\n'
            for k, v in self.Param.items():
                ret += str(k) + ' ' + str(v).strip('[] ') + '\n' 
            ret += '</' + self.Name + '>\n'
            return ret
