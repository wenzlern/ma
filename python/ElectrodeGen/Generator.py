#!/usr/bin/env python
"""This is the main file to call to generate a range of Electrodes"""

__author__     = "Nils Wenzler"
__copyright__  = "Nils Wenzler"
__maintainer__ = "Nils Wenzler"
__email__      = "wenzlern@ethz.ch"

# Imports
from PlateletClass import *
from ElectrodeGenerator import *
from ElectrodeClass import *
import argparse
from multiprocessing import Pool
import os
import matplotlib.pyplot as plt
import numpy as np

from mpl_toolkits.mplot3d import axes3d, Axes3D

# Functions 
def GetNrOfPlatelets(size, plateletdim, porosity):
    # Calculate the volume of the electrode and the platelet
    ElectrodeVolume = np.prod(size)
    PlateletVolume  = 0.5*np.prod(np.asarray(plateletdim))

    # Find the number of platelets we want
    NrOfPlatelets = int(np.round((1-porosity)*(ElectrodeVolume/PlateletVolume)))

    print 'Electrode Volume [um^3]: ', ElectrodeVolume
    print 'Platelet Volume [um^3]: ', PlateletVolume
    print 'NrOfPlatelets: ', NrOfPlatelets

    return NrOfPlatelets


def GetNewPlatelets(size, plateletdim, voxelsize, porosity, disparity,
                    NrOfPlatelets, color=None, align=None):
    Platelets = [GetRandomPlatelet(1, size, plateletdim, disparity=disparity, align=align)]
    for i in range(2,NrOfPlatelets+1):
        while True:
            checked = 0
            # Generate a random platelet
            temp = GetRandomPlatelet(i, size, plateletdim, disparity=disparity, color=color, align=align)
            for s in Platelets:
                if InterferenceDistance(s, temp, plateletdim):
                    break 
                else:
                    checked += 1
            if checked == len(Platelets):
                Platelets.append(temp)
                break

    print 'Generated  ', NrOfPlatelets, ' Platelets'
    return Platelets

def NonOverlapping(tup):
    TestedPlatelets, NewPlatelet = tup
    # Artifically set the new one to modified
    NewPlatelet.Modified = True
    while NewPlatelet.Modified == True:
        # Reset flag before commiting it to the test (which can alter the object)
        NewPlatelet.Modified=False

        for s in TestedPlatelets:
            if InterferenceRand(s, NewPlatelet, args.platelet, 10):
            #if InterferenceDistance(s, NewPlatelet, args.platelet):
                return

    # If we got until here we succeded
    return NewPlatelet


def AlignNonOverlapping(tup):
    TestedPlatelets, NewPlatelet = tup
    # Artifically set the new one to modified
    NewPlatelet.Modified = True
    while NewPlatelet.Modified == True:
        # Reset flag before commiting it to the test (which can alter the object)
        NewPlatelet.Modified=False

        for s in TestedPlatelets:
            if InterferenceAlign(s, NewPlatelet, args.platelet, 100):
                return

    # If we got until here we succeded
    return NewPlatelet


# Script
## Parse command line arguments
parser = argparse.ArgumentParser()
# Command line arguments
parser.add_argument('-o', '--output', nargs='?', type=str,
                    default = 'Electrode',
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
                    default=0.0, help='The disparity of particle size [0,1]')

parser.add_argument('-b', '--blocksize', nargs='?', type=float,
                    default=50, help='The blocksize to generate in single process mode')

parser.add_argument('-c', '--color', action='store_true',
                    help='Assign a random color to each platelet for illustration')

parser.add_argument('-a', '--align', action='store_true',
                    help='Align the platelet with the large face normal to x')

parser.add_argument('--stats', action='store_true',
                    help='Display statistics about the electrode')


args = parser.parse_args()


# Convert the dimensions into um from voxelcount
size = np.asarray(args.size)*args.voxel

NrOfPlatelets = GetNrOfPlatelets(size, args.platelet, args.porosity)

# Get a bunch of noninterfering particles to start with 
# We need to 'define' global variables here
TestedPlatelets = GetNewPlatelets(size, args.platelet, args.voxel, 
                                  args.porosity, args.disparity, args.blocksize, align=args.align)

# Get a worker pool
pool = Pool()

i = 1
# Now fill in the remaining particles
while len(TestedPlatelets) != NrOfPlatelets:
    print '############## Iteration ', i, ' ##############'
    # Remaining platelets to generate
    Remaining = min(NrOfPlatelets-len(TestedPlatelets), args.blocksize)
    # Generate platelets to fill in
    UntestedPlatelets = GetNewPlatelets(size, args.platelet, args.voxel, 
                            args.porosity, args.disparity, Remaining, color=args.color, align=args.align)
    if args.align:
        GoodPlatelets = pool.map(AlignNonOverlapping, [(TestedPlatelets, x) for x in UntestedPlatelets])
    else:
        GoodPlatelets = pool.map(NonOverlapping, [(TestedPlatelets, x) for x in UntestedPlatelets])
    print 'Found ', len(filter(None,GoodPlatelets)), ' non-interfering Platelets'
    # Append them to our list
    TestedPlatelets.extend(filter(None,GoodPlatelets))

    print len(TestedPlatelets), ' of ', NrOfPlatelets, ' Platelets done'

    # Just our Iteration counter
    i += 1

pool.close()

####### Debug
## Remaining platelets to generate
#Remaining = min(NrOfPlatelets-len(TestedPlatelets), args.blocksize)
## Generate platelets to fill in
#UntestedPlatelets = GetNewPlatelets(size, args.platelet, args.voxel, 
#                        args.porosity, args.disparity, Remaining, color=args.color)
#GoodPlatelets = pool.map(AlignNonOverlapping, UntestedPlatelets)
## Append them to our list
#TestedPlatelets.extend(filter(None,GoodPlatelets))
#
#print len(TestedPlatelets), ' of ', NrOfPlatelets, ' Platelets done'
#
#pool.close()
####### Debug


# Get a new electrode object
Anode = Electrode(args.size, args.voxel)

# We need to sanitize the numbers of the objects
for i in range(len(TestedPlatelets)):
    TestedPlatelets[i].SetNr(i+1)

Anode.Objects = TestedPlatelets
Anode.CreateMacro(args.output)
os.system('/usr/scratch/vanadium/lne-software/geodict/GeoDict2013/geodict2013 ' + args.output + '.gmc')

if args.stats:
    # Now we also want to display some stats
    # Calculate the distribution of directions
    # Our normal direction is the 100 (x-dir)
#    DotProducts = []
#    for s in Anode.Objects:
#        DotProducts.append(np.dot(s.Param['Axis1'], [0,1,0]))
#
#    Angles = 57.276*np.arccos(DotProducts)
#
#    print 'Mean of dot Product with 100: ', np.mean(DotProducts)
#    print 'Mean of angle with 100: ', np.mean(Angles)
#    Distribution = np.histogram(Angles, bins=np.linspace(0,180,num=21, endpoint=True))
#    plt.plot(np.linspace(0,180,20), Distribution[0].astype(float)/float(len(Distribution[0])))
    
#    
    x = []
    y = []
    z = []
    for s in Anode.Objects:
        x.append(s.Param['Axis1'][0])
        y.append(s.Param['Axis1'][1])
        z.append(s.Param['Axis1'][2])

#    fig = plt.figure()
#    ax = Axes3D(fig)
#    ax.set_aspect("equal")
#    
#    ax.contour(Axis[:][0], Axis[:][1], Axis[:][2])#, 16, extend3d=True)

    

    histx = np.histogram(x, bins=np.linspace(-1,1,num=21, endpoint=True))
    histy = np.histogram(y, bins=np.linspace(-1,1,num=21, endpoint=True))
    histz = np.histogram(z, bins=np.linspace(-1,1,num=21, endpoint=True))
    plt.plot(np.linspace(-1,1,20), histx[0].astype(float)/float(len(x)), label='x')
    plt.plot(np.linspace(-1,1,20), histy[0].astype(float)/float(len(y)), label='y')
    plt.plot(np.linspace(-1,1,20), histz[0].astype(float)/float(len(z)), label='z')
    plt.axis([-1,1,0,1])


#    Distribution = np.histogram(Axis, bins=np.linspace(-1,1,num=31, endpoint=True))
#    plt.plot(np.linspace(-1,1,20), Distribution[0].astype(float)/float(len(Distribution[0])))
    plt.legend()
    plt.show()
