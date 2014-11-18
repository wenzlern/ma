"""Functions to map vtk points to cells"""

__author__     = "Nils Wenzler"
__copyright__  = "Nils Wenzler"
__maintainer__ = "Nils Wenzler"
__email__      = "wenzlern@ethz.ch"

import numpy as np

def SanitizeCellArray(cells):
  """Sanitize the cell list"""
  cell_list = []
  for i in range(0, len(cells)/9, 9):
    cell_list.append(sorted(cells[i+1:i+9]))

  return cell_list
 
class CellIndices:
  """Class providing a CellIndices matcher"""

  def __init__(self,points, cells, PointGridSpacing):
    self.PointList = points//PointGridSpacing
    self.CellPoints = np.zeros(8)
    self.CellList = SanitizeCellArray(cells)

  def FindPointIndex(self,x, y, z):
    """Macro to find index in array"""
    return self.PointList.tolist().index([x,y,z])
  
 
  def FindCellIndex(self, x, y, z, ):
    """Find the index of the cell provided a location x, y, z"""
  
    #These are the eight points making our cell
    self.CellPoints[0] = self.FindPointIndex(x  , y  , z  )
    self.CellPoints[1] = self.FindPointIndex(x+1, y  , z  )
    self.CellPoints[2] = self.FindPointIndex(x+1, y+1, z  )
    self.CellPoints[3] = self.FindPointIndex(x  , y+1, z  )
    self.CellPoints[4] = self.FindPointIndex(x  , y  , z+1)
    self.CellPoints[5] = self.FindPointIndex(x+1, y  , z+1)
    self.CellPoints[6] = self.FindPointIndex(x+1, y+1, z+1)
    self.CellPoints[7] = self.FindPointIndex(x  , y+1, z+1)

    CellIndex = self.CellList.index(sorted(self.CellPoints))
  
    return CellIndex
  
  def FindCellIndices(self, xrng, yrng, zrng):
    """Top level function returning a list for a given x, y, z range"""
    Indices = []
    
    #Loop over the range
    for z in zrng:
      for y in yrng:
        for x in xrng:
          Indices.append(self.FindCellIndex(x, y, z))
          print Indices[-1]
  
