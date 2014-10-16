"""
@summary: This module contains a base class for defining a group of compression
             algorithms that use a two-stage approach by first using a tree
             compression technique until a threshold is met, then switch to a 
             run-length encoding method.
@author: CJ Grady
@version: 2.0
@status: beta

@license: gpl2
@copyright: Copyright (C) 2014, University of Kansas Center for Research

          Lifemapper Project, lifemapper [at] ku [dot] edu, 
          Biodiversity Institute,
          1345 Jayhawk Boulevard, Lawrence, Kansas, 66045, USA
   
          This program is free software; you can redistribute it and/or modify 
          it under the terms of the GNU General Public License as published by 
          the Free Software Foundation; either version 2 of the License, or (at 
          your option) any later version.
  
          This program is distributed in the hope that it will be useful, but 
          WITHOUT ANY WARRANTY; without even the implied warranty of 
          MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
          General Public License for more details.
  
          You should have received a copy of the GNU General Public License 
          along with this program; if not, write to the Free Software 
          Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
          02110-1301, USA.
"""

from matrix.matrix import _CompressedGrid

# .............................................................................
class _ComboCompressionBaseClass(_CompressedGrid):
   """
   @summary: This is the base class for two-stage compression algorithms
   """
   # .............................
   def __init__(self, rleMethod, orderThreshold=8, grid=None):
      """
      @summary: Base class constructor for combination (two-stage) compression
      @param rleMethod: This is the run-length encoding compression method to 
                           use once the threshold is reached
      @param orderThreshold: (optional: default=8), once there are fewer than
                                two to this number elements (2^orderThreshold)
                                in the grid section, switch to the run-length
                                encoding method.
      @param grid: (optional) A Grid to compress 
      """
      self.rleMethod = rleMethod
      self.threshold = 2**orderThreshold
      if grid is not None:
         self.compress(grid)
      else:
         self._initialize()
   
   # .............................
   def _initialize(self):
      """
      @summary: Use this method in sub-classes to initialize any variables
      """
      pass
   