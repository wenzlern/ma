"""Electrode Storage Class
This class stores all the information we get about the electrode
out of a BEST simulation. It combines all the output files into one
object with the data arrays spanning over time. Later efficient binary
store and load operations are given"""

__author__     = "Nils Wenzler"
__copyright__  = "Nils Wenzler"
__maintainer__ = "Nils Wenzler"
__email__      = "wenzlern@ethz.ch"

import ctypes
from enum import Enum

class Materials(Enum):
  """Provide material identifier by name"""
  #range(1,3) generates a list [1,2], therefor + 1
  #Ranges can be found in BEST parametrization
  Active = range(100, 110 + 1)
  Electrolyte = range(200, 210 + 1)


class Cell:
  """Stores the individual Cell information"""
  def __init__(self, MaterialID, Concentration, Potential, OverPotential,
               IndividualOCV, ParticleFlux, CurrentDensity):
    # Scalars
    self.MaterialID    = MaterialID    
    self.Concentration = Concentration 
    self.Potential     = Potential     
    self.OverPotential = OverPotential 
    self.IndividualOCV = IndividualOCV 
    # Vectors
    self.ParticleFlux   = ParticleFlux  
    self.CurrentDensity = CurrentDensity


class Electrode:
  """Store the geometry and access the cells"""
  #Create an object that contains pointers to cells
  #in format [x,y,z]
  def __init__(self, xsize, ysize, zsize, timestep = 30):
    """Initizialize Electrode object with size"""
    #Creates an list of lists of lists of None references
    self.Cells = [[[None for z in range(zsize)] 
                 for y in range(ysize)] 
                 for x in range(xsize)]
    self.TimeStep = timestep


  def SetCell(self, x, y, z, MaterialID, Concentration, Potential, 
              OverPotential, IndividualOCV, ParticleFlux, CurrentDensity):
      """Set the cell at position x,y,z to the provided values"""
      self.Cells[x][y][z] = Cell(MaterialID, Concentration, Potential, 
                                 OverPotential, IndividualOCV, ParticleFlux, 
                                 CurrentDensity)
       

  #The timestep from dataset to dataset

  #The size of the electrode
