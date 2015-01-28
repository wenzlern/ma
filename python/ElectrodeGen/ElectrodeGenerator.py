#!/usr/bin/env python
"""This is the main file to call to generate a range of Electrodes"""

__author__     = "Nils Wenzler"
__copyright__  = "Nils Wenzler"
__maintainer__ = "Nils Wenzler"
__email__      = "wenzlern@ethz.ch"

# Imports
from PlateletClass import *
import random
import sys

def statusbar(progress, total):  
    print('\r[{0:10}]{1:>2}%'.format('#' * int(progress * 10 /total), progress))

def GetRandomPlatelet(number, electrodedim, plateletdim, disparity=None, color=None):
    position = [random.uniform(0, electrodedim[0]), random.uniform(0, electrodedim[1]),
                random.uniform(0, electrodedim[2])]
    direction = [random.uniform(-1, 1), random.uniform(-1, 1),
                random.uniform(-1, 1)]
    angle = random.randint(0,360)

    return Platelet(number, position, direction, angle,
                    longdiam = plateletdim[0], shortdiam = plateletdim[1],
                    thickness = plateletdim[2], sizedisp = disparity, color=color)

def Interference(object1, object2, objectdim):
    dist = np.linalg.norm(object1.GetPos()-object2.GetPos())
    object1Ax3 = Rotate3DVector(object1.GetAx2(), object1.GetAx1(), 90)
    object2Ax3 = Rotate3DVector(object2.GetAx2(), object2.GetAx1(), 90)
    if dist < (objectdim[0]):
        if dist < (np.dot(object1.GetAx2(), object2.GetAx2())*objectdim[0]):
            return True 
        elif dist < (np.dot(object1Ax3, object2Ax3)*objectdim[1]):
            return True
        else:
            return False
    else:
        return False

def InterferenceDistance(object1, object2, objectdim):
    dist = np.linalg.norm(object1.GetPos()-object2.GetPos())
    if dist < (objectdim[0]):
        return True 
    else:
        return False

def InterferenceAlign(object1, object2, objectdim, NrOfPos):
    if NrOfPos == 0:
        return True
    else:
        dist = np.linalg.norm(object1.GetPos()-object2.GetPos())
        object1Ax3 = Rotate3DVector(object1.GetAx2(), object1.GetAx1(), 90)
        object2Ax3 = Rotate3DVector(object2.GetAx2(), object2.GetAx1(), 90)
        if dist < (objectdim[0]):
            if (dist < (np.dot(object1.GetAx2(), object2.GetAx2())*objectdim[0])) or (dist < (np.dot(object1Ax3, object2Ax3)*objectdim[1])) or (dist < (np.dot(-object1.GetAx2(), object2.GetAx2())*objectdim[0])) or (dist < (np.dot(-object1Ax3, object2Ax3)*objectdim[1])) or (dist < (np.dot(object1.GetAx2(), -object2.GetAx2())*objectdim[0])) or (dist < (np.dot(object1Ax3, -object2Ax3)*objectdim[1])) or (dist < (np.dot(-object1.GetAx2(), -object2.GetAx2())*objectdim[0])) or (dist < (np.dot(-object1Ax3, -object2Ax3)*objectdim[1])):
                ax1, ax2 = GetAxis(object1.GetAx1(), NrOfPos)
                #object2.SetAx1([ax1[0]*random.uniform(0,1), 
                #                  ax1[1]*random.uniform(0,1), 
                #                  ax1[1]*random.uniform(0,1)])

                #object2.SetAx2([ax2[0]*random.uniform(0,1), 
                #                  ax2[1]*random.uniform(0,1), 
                #                  ax2[2]*random.uniform(0,1)])

                object2.SetAx1([random.uniform(-1,1), 
                                random.uniform(-1,1), 
                                random.uniform(-1,1)])

                object2.SetAx2([random.uniform(-1,1), 
                                random.uniform(-1,1), 
                                random.uniform(-1,1)])

                return InterferenceAlign(object1, object2, objectdim, NrOfPos - 1)
            else:
                return False
        else:
            return False


# Generators
def RandomOverlapping(size, plateletdim, voxelsize, 
                      porosity, disparity):

    # Calculate the volume of the electrode and the platelet
    ElectrodeVolume = np.prod(size)
    PlateletVolume  = 0.5*np.prod(np.asarray(plateletdim))

    # Find the number of platelets we want
    NrOfPlatelets = int(np.round((1-porosity)*(ElectrodeVolume/PlateletVolume)))


    print 'Electrode Volume [um^3]: ', ElectrodeVolume
    print 'Platelet Volume [um^3]: ', PlateletVolume
    print 'NrOfPlatelets: ', NrOfPlatelets

    Platelets = []
    for i in range(NrOfPlatelets):
        Platelets.append(GetRandomPlatelet(i+1, size, plateletdim, disparity=disparity))

    return Platelets

def RandomNonOverlapping(plateletdim, voxelsize, 
                      porosity, disparity):

    # Calculate the volume of the electrode and the platelet
    ElectrodeVolume = np.prod(size)
    PlateletVolume  = 0.5*np.prod(np.asarray(plateletdim))

    # Find the number of platelets we want
    NrOfPlatelets = int(np.round((1-porosity)*(ElectrodeVolume/PlateletVolume)))

    print 'Electrode Volume [um^3]: ', ElectrodeVolume
    print 'Platelet Volume [um^3]: ', PlateletVolume
    print 'NrOfPlatelets: ', NrOfPlatelets

    Platelets = [GetRandomPlatelet(1, size, plateletdim, disparity=disparity)]
    for i in range(2,NrOfPlatelets+1):
        while True:
            checked = 0
            # Generate a random platelet
            temp = GetRandomPlatelet(i, size, plateletdim, disparity=disparity)
            for s in Platelets:
                if Interference(s, temp, plateletdim):
                    break 
                else:
                    checked += 1
            if checked == len(Platelets):
                Platelets.append(temp)
                break

        print 'Platetelet ', i, 'of ', NrOfPlatelets

    return Platelets

def AlignNonOverlapping(size, plateletdim, voxelsize, 
                         porosity, disparity):

    # Calculate the volume of the electrode and the platelet
    ElectrodeVolume = np.prod(size)
    PlateletVolume  = 0.5*np.prod(np.asarray(plateletdim))

    # Find the number of platelets we want
    NrOfPlatelets = int(np.round((1-porosity)*(ElectrodeVolume/PlateletVolume)))

    print 'Electrode Volume [um^3]: ', ElectrodeVolume
    print 'Platelet Volume [um^3]: ', PlateletVolume
    print 'NrOfPlatelets: ', NrOfPlatelets

    Platelets = [GetRandomPlatelet(1, size, plateletdim, disparity=disparity)]
    for i in range(2,NrOfPlatelets+1):
        while True:
            checked = 0
            # Generate a random platelet
            temp = GetRandomPlatelet(i, size, plateletdim, disparity=disparity)
            for s in Platelets:
                if InterferenceAlign(s, temp, plateletdim, 180):
                    break 
                else:
                    checked += 1
            if checked == len(Platelets):
                Platelets.append(temp)
                break

        print 'Platetelet ', i, 'of ', NrOfPlatelets

    return Platelets

# If called directly
if __name__ == "__main__":
    from ElectrodeClass import *
    import argparse

    ## Parse command line arguments
    parser = argparse.ArgumentParser()
    # Command line arguments
    parser.add_argument('-o', '--output', nargs='?', type=str,
                        default = 'Electrode.gps',
                        help='The name/path of the output file')
    
    parser.add_argument('-s', '--size', nargs='*', type=int,
                        help='The size of the electrode in voxels {x y z}')
    
    parser.add_argument('-l', '--platelet', nargs='*', type=float,
                        help='The size of the platelet in um {ld, sd, t}\n' + 
                        'ld := long diameter, sd := short diameter, t := thickness')
    
    parser.add_argument('-v', '--voxel', nargs='?', type=float,
                        default=5e-7, help='The size of one Voxel')
    
    parser.add_argument('-p', '--porosity', nargs='?', type=float,
                        default=0.3, help='The porosity of the Electrode [0,1]')
    
    parser.add_argument('-d', '--disparity', nargs='?', type=float,
                        default=0.1, help='The disparity of particle size [0,1]')
    
    args = parser.parse_args()
    
    
    # Get a new electrode object
    Anode = Electrode(args.size, args.voxel)

    Anode.Objects = AlignNonOverlapping(np.asarray(args.size)*args.voxel, args.platelet, args.voxel, 
                                         args.porosity, args.disparity)

    dir = os.getcwd()

    Anode.SaveToFile(dir + args.output)
