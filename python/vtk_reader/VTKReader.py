"""VTK files reader functions
This class aims to simplyfie the reading in of vtk files generated 
by the BEST tools."""

__author__     = "Nils Wenzler"
__copyright__  = "Nils Wenzler"
__maintainer__ = "Nils Wenzler"
__email__      = "wenzlern@ethz.ch"

##imports
import sys
import argparse
import glob
import os
import numpy as np
import scipy.io as spio
import xml.etree.ElementTree as et

from vtk import *
from vtk.util.numpy_support import vtk_to_numpy

from PointsToCell import *


##Some general functions
#Run a Match
def RunMatch(points, cells, PointGridSpacing, xrng, yrng, zrng):
    MatchIndices = CellIndices(points, cells, PointGridSpacing)
    return MatchIndices.FindCellIndices(xrng, yrng, zrng)

#Reorder Scalars
def VTKScalarToXYZ(Array, xsize, ysize, zsize):
    NumpyArray = vtk_to_numpy(Array)
    return NumpyArray.reshape((xsize,ysize,zsize), order='F')

#Reorder arrays
def VTKVectorToXYZ(Array, xsize, ysize, zsize):
    NumpyArray = vtk_to_numpy(Array)
    return NumpyArray.reshape((xsize,ysize,zsize,3), order='F')

#Sorting function
def FileNumber(s):
    return int(str.split(str.split(s,'_')[-1],'.')[0])



## Parse command line arguments
parser = argparse.ArgumentParser()
# Command line arguments
parser.add_argument('foldername', nargs='?', type=str,  
                    help='The Simulation output folder name')

parser.add_argument('-d', '--destination', nargs='+', type=str,
                    help='The destination for the array files to be saved.'
                          + 'Default is the simulation output folder')

parser.add_argument('-c', '--check', action='store_true',
                    default=False, help='Enable the checking of the'
                    +'geometry. Usually not necessary, but do once on'
                    +'new geometries')

parser.add_argument('-g', '--gridspacing', nargs='?', type=float,
                    default=6.5e-5, help='Provide the grid spacing,' 
                    +'default is 6.5e-5')

parser.add_argument('-s', '--dimension', nargs='*', type=int,
                      help='The size of the electrode in voxels {x y z}')

parser.add_argument('-t', '--timestep', nargs='?', type=int, default=30,
                      help='The timestep of the vtk writeout, default is 30')

parser.add_argument('-m', '--matlab', action='store_true', default=False, 
                      help='Store arrays in matlab format (.mat),'
                            +'default is numpy (.npy)')

parser.add_argument('-z', '--compress', action='store_true', default=False, 
                      help='Store in compresses .npz file, default is false')

args = parser.parse_args()

##Sanitize input
#Default dest is input folder
if not args.destination:
    args.destination = args.foldername

#For convenience store sizes to vars
xsize = args.dimension[0]
ysize = args.dimension[1]
zsize = args.dimension[2]

#Get the files to process and sort them by number
Files = []
for i in glob.iglob(os.path.join(args.foldername + '/output/', '*.vtk')):
    Files.append(i)

Files.sort(key=FileNumber)

##Main
#Empty lists for easy appending of all vtk data
Concentration      = [] 
Potential          = [] 
Overpotential      = [] 
MaterialIdentifier = []
IndividualOCV      = [] 
ParticleFlux       = [] 
CurrentDensity     = []
FileIndex          = []

#Get reader for UnstructuredGrid
reader = vtkUnstructuredGridReader()

#Go through all files
for filename in Files:
    print 'Reading ' + filename
    reader.SetFileName(filename)
    reader.ReadAllScalarsOn()
    reader.ReadAllVectorsOn()
    reader.ReadAllFieldsOn()
    reader.Update()
    ugrid = reader.GetOutput()
    
    #Only perform this check if specified on one file
    if FileNumber(filename) == 0 and args.check == True:
        points = ugrid.GetPoints().GetData()
        cells = ugrid.GetCells().GetData()
        
        points = vtk_to_numpy(points)
        cells  = vtk_to_numpy(cells)
        
        indices = RunMatch(points, cells, args.gridspacing, 
                           range(xsize), range(ysize), range(zsize))
        
        #If nothing changes when sorting, it is ok
        if indices.sort() == indices:
            print "Geometry passed check, Format as expected"
        else:
            print "Geometry did not pass, different Format as expected"
            print "You might want to adapt the script to reflect the changes"
            sys.exit()
         
    #Read out Scalars
    Concentration.append(VTKScalarToXYZ(ugrid.GetCellData().GetScalars("concentration"),
                                     xsize, ysize, zsize))
    Potential.append(VTKScalarToXYZ(ugrid.GetCellData().GetScalars("potential"),
                                     xsize, ysize, zsize))
    Overpotential.append(VTKScalarToXYZ(ugrid.GetCellData().GetScalars("overpotential"),
                                     xsize, ysize, zsize))
    MaterialIdentifier.append(VTKScalarToXYZ(ugrid.GetCellData().GetScalars("materialIdentifier"),
                                     xsize, ysize, zsize))
    IndividualOCV.append(VTKScalarToXYZ(ugrid.GetCellData().GetScalars("individualOCV"),
                                     xsize, ysize, zsize))
    
    #Read out Vectors
    ParticleFlux.append(VTKVectorToXYZ(ugrid.GetCellData().GetVectors("particle_flux"),
                         xsize, ysize, zsize))
    CurrentDensity.append(VTKVectorToXYZ(ugrid.GetCellData().GetVectors("current_density"),
                         xsize, ysize, zsize))
    



#Convert lists to arrays
FileIndex = np.asarray(FileIndex)
Concentration = np.asarray(Concentration)
Potential = np.asarray(Potential)
Overpotential = np.asarray(Overpotential)
MaterialIdentifier = np.asarray(MaterialIdentifier)
IndividualOCV = np.asarray(IndividualOCV)
ParticleFlux = np.asarray(ParticleFlux)
CurrentDensity = np.asarray(CurrentDensity)


#We also would like to have the parameters and the paramters_application
Parameters = et.parse(args.foldername + '/configuration/' + 'parameters.xml').getroot()
ParamApplication = et.parse(args.foldername + '/configuration/' + 'parameters_application.xml').getroot()


##Store the arrays into npy, npz or mat files
#Extract File name
Filename = str.split(args.foldername,'/')[-1]

if args.matlab:
    spio.savemat(args.destination + '/' + Filename + '.mat',
               {'Concentration':Concentration,
                'Potential':Potential,
                'Overpotential':Overpotential, 
                'MaterialIdentifier':MaterialIdentifier,
                'IndividualOCV':IndividualOCV,
                'ParticleFlux':ParticleFlux,
                'CurrentDensity':CurrentDensity})
    print 'Storing to ' + args.destination + '/' + Filename + '.mat'

elif args.compress:
    np.savez_compressed(args.destination + '/' + Filename, Concentration=Concentration,
                      Potential=Potential, OverPotential=Overpotential,
                      MaterialIdentifier=MaterialIdentifier, IndividualOCV=IndividualOCV,
                      ParticleFlux=ParticleFlux, CurrentDensity=CurrentDensity, Parameters=Parameters,
                      ParamApplication=ParamApplication)
    print 'Storing to ' + args.destination + '/' + Filename + 'npz'

else:
    np.savez(args.destination + '/' + Filename, Concentration=Concentration,
           Potential=Potential, OverPotential=Overpotential,
           MaterialIdentifier=MaterialIdentifier, IndividualOCV=IndividualOCV,
           ParticleFlux=ParticleFlux, CurrentDensity=CurrentDensity, Parameters=Parameters,
           ParamApplication=ParamApplication)

    print 'Storing to ' + args.destination + '/' + Filename + '.npy'

