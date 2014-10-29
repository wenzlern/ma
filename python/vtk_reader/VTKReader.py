"""VTK files reader functions
This class aims to simplyfie the reading in of vtk files generated 
by the BEST tools."""

__author__     = "Nils Wenzler"
__copyright__  = "Nils Wenzler"
__maintainer__ = "Nils Wenzler"
__email__      = "wenzlern@ethz.ch"

#imports
import sys
import argparse
from vtk import *
import glob
import os
from vtk.util.numpy_support import vtk_to_numpy
import numpy as np
from PointsToCell import *
import scipy.io as spio
import tables

def RunMatch(points, cells, PointGridSpacing, xrng, yrng, zrng):
  MatchIndices = CellIndices(points, cells, PointGridSpacing)
  return MatchIndices.FindCellIndices(xrng, yrng, zrng)

def VTKToXYZ(Array, xsize, ysize, zsize):
  NumpyArray = vtk_to_numpy(Array)
#  XYArray = []
#  XYZArray = []
#  #Reassemble in 2 steps to a 3D list
#  for y in range(ysize):
#    XYArray.append([y*xsize:((y+1)*xsize)-1])
#  for z in range(zsize):
#    XYZArray.append(XYArray[z*ysize:((z+1)*ysize)-1])
  return NumpyArray.reshape((xsize,ysize,zsize), order='F')



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

#Sanitize input
if not args.destination:
  args.destination = args.foldername

#For convenience store sizes
xsize = args.dimension[0]
ysize = args.dimension[1]
zsize = args.dimension[2]

#Get the number of files to process
#NrFiles = len([item for item in os.listdir(args.foldername + '/output') if os.path.isfile(os.path.join(args.foldername + '/output', item))])
#print str(NrFiles) + ' files to process'

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

for filename in glob.iglob(os.path.join(args.foldername + '/output/', '*.vtk')):
  print 'Reading ' + filename
  reader.SetFileName(filename)
  reader.ReadAllScalarsOn()
  reader.ReadAllVectorsOn()
  reader.ReadAllFieldsOn()
  reader.Update()
  ugrid = reader.GetOutput()

  if filename == 'output_default_0.vtk' and args.check == True:
    points = ugrid.GetPoints().GetData()
    cells = ugrid.GetCells().GetData()
    
    points = vtk_to_numpy(points)
    cells  = vtk_to_numpy(cells)
    
    indices = RunMatch(points, cells, args.gridspacing, 
                       range(xsize), range(ysize), range(zsize))
    
    if indices.sort() == indices:
      print "Geometry passed check, Format as expected"
    else:
      print "Geometry did not pass, different Format as expected"
      print "You might want to adapt the script to reflect the changes"
      sys.exit()
     
  #Get Index to write to from filename
  FileIndex.append(int(str.split(str.split(filename,'.')[0],'_')[-1]))

  #Read out Scalars
  Concentration.append(VTKToXYZ(ugrid.GetCellData().GetScalars("concentration"),
                                   xsize, ysize, zsize))
  Potential.append(VTKToXYZ(ugrid.GetCellData().GetScalars("potential"),
                                   xsize, ysize, zsize))
  Overpotential.append(VTKToXYZ(ugrid.GetCellData().GetScalars("overpotential"),
                                   xsize, ysize, zsize))
  MaterialIdentifier.append(VTKToXYZ(ugrid.GetCellData().GetScalars("materialIdentifier"),
                                   xsize, ysize, zsize))
  IndividualOCV.append(VTKToXYZ(ugrid.GetCellData().GetScalars("individualOCV"),
                                   xsize, ysize, zsize))

  #Read out Vectors
#  ParticleFlux   = VTKToXYZ(ugrid.GetCellData().GetVectors("particle_flux"),
#                       xsize, ysize, zsize)
#  CurrentDensity = VTKToXYZ(ugrid.GetCellData().GetVectors("current_density"),
#                       xsize, ysize, zsize)


#Convert lists to arrays and Sort arrays according to FileIndex
FileIndex = np.asarray(FileIndex)
Concentration = np.asarray(Concentration)
Potential = np.asarray(Potential)
Overpotential = np.asarray(Overpotential)
MaterialIdentifier = np.asarray(MaterialIdentifier)
IndividualOCV = np.asarray(IndividualOCV)
ParticleFlux = np.asarray(ParticleFlux)
CurrentDensity = np.asarray(CurrentDensity)

#Dictionary of arrays to store
#arrays = {'Concentration'     : Concentration,
#          'Potential'         : Potential,
#          'Overpotential'     : Overpotential, 
#          'MaterialIdentifier': MaterialIdentifier,
#          'IndividualOCV'     : IndividualOCV} 
#          'ParticleFlux'      : ParticleFlux,
#          'CurrentDensity'    : CurrentDensity}

#Extract File name
Filename = str.split(args.foldername,'/')[-1]

#Store the output files, we use Hfd5, since both matlab and python can use it
#FileHandler = h5py.File(args.destination + '/' + Filename + '.h5', 'w')
#FileHandler = tables.openFile(args.destination + '/' + Filename + '.h5', 'w')
#root = FileHandler.root
#
#FileHandler.createArray(root, 'Concentration', Concentration)
#FileHandler.createArray(root, 'Potential', Potential)
#FileHandler.createArray(root, 'OverPotential', Overpotential)
#FileHandler.createArray(root, 'MaterialIdentifier', MaterialIdentifier)
#FileHandler.createArray(root, 'IndividualOCV', IndividualOCV)


#Create datasets
#FileHandler.create_dataset('Concentration', data=Concentration)
#FileHandler.create_dataset('Potential', data=Potential)
#FileHandler.create_dataset('OverPotential', data=Overpotential)
#FileHandler.create_dataset('MaterialIdentifier', data=MaterialIdentifier)
#FileHandler.create_dataset('IndividualOCV', data=IndividualOCV)

#FileHandler.close()

if args.matlab:
  spio.savemat(args.destination + '/' + Filename + '.mat',
               {'Concentration':Concentration,
                'Potential':Potential,
                'Overpotential':Overpotential, 
                'MaterialIdentifier':MaterialIdentifier,
                'IndividualOCV':IndividualOCV})

elif args.compress:
  np.savez_compressed(args.destination + '/' + Filename, Concentration=Concentration,
                      Potential=Potential, OverPotential=Overpotential,
                      MaterialIdentifier=MaterialIdentifier, IndividualOCV=IndividualOCV,
                      ParticleFlux=ParticleFlux, CurrentDensity=CurrentDensity)
  print 'Storing to ' + args.destination + '/' + Filename + 'npz'

else:
  np.savez(args.destination + '/' + Filename, Concentration=Concentration,
           Potential=Potential, OverPotential=Overpotential,
           MaterialIdentifier=MaterialIdentifier, IndividualOCV=IndividualOCV,
           ParticleFlux=ParticleFlux, CurrentDensity=CurrentDensity)

  print 'Storing to ' + args.destination + '/' + Filename + '.npy'
  
#  #Now lets sort them into our electrode object
#  Electrode.append(args.dimension[0], args.dimension[1], args.dimension[2],
#                   args.timestep)
#  for z in range(args.dimension[2]):
#     for y in range(args.dimension[1]):
#        for x in range(args.dimension[0]:
#          #Calculate the index: we step through x then y then z
#          i = z*(
#          Electrode.SetCell(x,y,z, 
  


