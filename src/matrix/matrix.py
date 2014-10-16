"""
@summary: Module contain matrix base classes
@author: CJ Grady
@version: 1.0
@status: alpha

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

class Grid(object):
   """
   @summary: Base class for Lifemapper grids.  This class can be used with 
                uncompressed grids.
   """
   # ...........................
   def __init__(self, griddedData=None):
      if griddedData is not None:
         self._initFromGrid(griddedData)
      else:
         self.ySize = None
         self.xSize = None
         self.data = []
         self.classes = set([])

   # ...........................
   def _initFromGrid(self, griddedData):
      self.ySize = len(griddedData)
      self.xSize = len(griddedData[0])
      self.data = griddedData
      self.findClasses()
      
   # ...........................
   def findClasses(self):
      """
      @summary: Finds all of the unique classes in the data
      """
      self.classes = set([])
      for row in self.data:
         for col in row:
            self.classes.add(col)
   
   # ...........................
   def query(self, x, y):
      return self.data[y][x]

   # ...........................
   def read(self, fn):
      self.data = []
      with open(fn) as f:
         for line in f:
            self.data.append([int(i) for i in line.split(' ')])

   # ...........................
   def write(self, fn):
      with open(fn, 'w') as f:
         for row in self.data:
            f.write('%s\n' % ' '.join([str(i) for i in row]))

# .............................................................................
class _CompressedGrid(Grid):
   # ...........................
   def __init__(self):
      raise Exception, "init must be implemented in sub class"   
   
   # ...........................
   def query(self, x, y):
      raise Exception, "Query must be implemented in sub class"
   
   # ...........................
   def read(self, fn):
      raise Exception, "Read must be implemented in sub class"   

   # ...........................
   def write(self, fn):
      raise Exception, "Write must be implemented in sub class"
      
