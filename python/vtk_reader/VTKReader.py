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
from ClassElectrode import *
import scipy

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
parser.add_argument('foldername', nargs='+', type=str,  
                    help='The Simulation output folder name')

parser.add_argument('-d', '--destination', nargs='?', type=str,
                    help='The destination for the array files to be saved.
                          Default is the simulation output folder')

parser.add_argument('-c', '--check', nargs='?', action='store_true',
                    default=False, help='Enable the checking of the 
                    geometry. Usually not necessary, but do once on 
                    new geometries')

parser.add_argument('-g', '--gridspacing', nargs='?', type=float,
                    default=6.5e-5, help='Provide the grid spacing, 
                    default is 6.5e-5')

parser.add_argument('-d', '--dimension', nargs='*', type=int,
                      help='The size of the electrode in voxels {x y z}')

parser.add_argument('-t', '--timestep', nargs='?', type=int, default=30
                      help='The timestep of the vtk writeout, default is 30')

parser.add_argument('-m', '--matlab', nargs='?', action='store_true', default=False 
                      help='Store arrays in matlab format (.mat),
                            default is numpy (.npy)')

parser.add_argument('-z', '--compress', nargs='?', action='store_true', default=False 
                      help='Store in compresses .npz file, default is false')

args = parser.parse_args()


#For convenience store sizes
xsize = args.dimension[0]
ysize = args.dimension[1]
zsize = args.dimension[2]

for filename in glob.iglob(os.path.join(folder, '/output/*.vtk')):
  reader.SetFileName(filename)
  reader.ReadAllScalarsOn()
  reader.ReadAllVectorsOn()
  reader.ReadAllFieldsOn()
  reader.Update()
  ugrid = reader.GetOutput()

  if filename == 'output_default_0.vtk' && args.check == True:
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
     

  #Read out Scalars
  Concentration      = VTKToXYZ(ugrid.GetCellData().GetScalars("concentration"),
                       xsize, ysize, zsize)
  Potential          = VTKToXYZ(ugrid.GetCellData().GetScalars("potential"),
                       xsize, ysize, zsize)
  Overpotential      = VTKToXYZ(ugrid.GetCellData().GetScalars("overpotential"),
                       xsize, ysize, zsize)
  MaterialIdentifier = VTKToXYZ(ugrid.GetCellData().GetScalars("MaterialIdentifier"),
                       xsize, ysize, zsize)
  IndividualOCV      = VTKToXYZ(ugrid.GetCellData().GetScalars("IndividualOCV"),
                       xsize, ysize, zsize)

  #Read out Vectors
  ParticleFlux   = VTKToXYZ(ugrid.GetCellData().GetScalars("particle_flux"),
                       xsize, ysize, zsize)
  CurrentDensity = VTKToXYZ(ugrid.GetCellData().GetScalars("current_density"),
                       xsize, ysize, zsize)

  #Store the output files
  if args.matlab:
    scipy.io.savemat(args.foldername + '/' + args.foldername +'.mat',
                     dict(Concentration=Concentration, Potential=Potential,
                          Overpotential=Overpotential, MaterialIdentifier=MaterialIdentifier,
                          IndividualOCV=IndividualOCV, ParticleFlux=ParticleFlux,
                          CurrentDensity=CurrentDensity))
  else if args.compressed:
    savez_compressed(args.foldername + '/' + args.foldername +'.npz',
                     Concentration=Concentration, Potential=Potential,
                     Overpotential=Overpotential, MaterialIdentifier=MaterialIdentifier,
                     IndividualOCV=IndividualOCV, ParticleFlux=ParticleFlux,
                     CurrentDensity=CurrentDensity))

  else:
    savez(args.foldername + '/' + args.foldername +'.npy',
          Concentration=Concentration, Potential=Potential,
          Overpotential=Overpotential, MaterialIdentifier=MaterialIdentifier,
          IndividualOCV=IndividualOCV, ParticleFlux=ParticleFlux,
          CurrentDensity=CurrentDensity))


#  #Now lets sort them into our electrode object
#  Electrode.append(args.dimension[0], args.dimension[1], args.dimension[2],
#                   args.timestep)
#  for z in range(args.dimension[2]):
#     for y in range(args.dimension[1]):
#        for x in range(args.dimension[0]:
#          #Calculate the index: we step through x then y then z
#          i = z*(
#          Electrode.SetCell(x,y,z, 
  


