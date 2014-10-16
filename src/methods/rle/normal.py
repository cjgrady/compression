"""
@summary: This module contains a class for compressing a grid using left to 
             right, top to bottom, run length encoding
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
import struct

from matrix.matrix import _CompressedGrid, Grid

# .............................................................................
class NormalRLECompressedGrid(_CompressedGrid):
   """
   @summary: This class compresses a grid using left to right run length 
                encoding
   """
   METHOD = 1
   VERSION = 2.0
   # ...........................
   def __init__(self, grid=None):
      if grid is not None:
         self.compress(grid)
      else:
         self.lyrs = []
   
   # ...........................
   def compress(self, grid):
      """
      @summary: Compresses a Grid into a NormalRLECompressedGrid
      @param grid: The Grid object to compress
      """
      self.xSize = grid.xSize
      self.ySize = grid.ySize
      
      lArray = []
      for row in grid.data:
         lArray.extend(row)
      
      self.cmpData = []
      
      # Compress by twos
      cur = (lArray[0], lArray[1])
      num = 1
      for i in xrange(2, len(lArray), 2):
         try:
            if cur == (lArray[i], lArray[i+1]):
               num += 1
            else:
               self.cmpData.append((cur, num))
               cur = (lArray[i], lArray[i+1])
               num = 1
         except IndexError:
            # This happens if there is an odd number of items
            self.cmpData.append((cur, num))
            cur = (lArray[-1])
            num = 1
      self.cmpData.append((cur, num))

   # ...........................
   def decompress(self):
      """
      @summary: Decompresses the compressed grid into a Grid object
      """
      data = []
      for y in xrange(self.ySize):
         row = []
         for x in xrange(self.xSize):
            row.append(0)
         data.append(row)
      
      lArray = []
      for cur, num in self.cmpData:
         lArray.extend(num*list(cur))
         
      x = 0
      y = 0
      for i in lArray:
         data[y][x] = i
         x += 1
         if x == self.xSize:
            x = 0
            y += 1
      
      return Grid(griddedData = data)
   
   # ...........................
   def query(self, x, y):
      """
      @summary: Queries the compressed grid to find the value at the specified 
                   coordinates
      @param x: The x (horizontal) coordinate, starts from the left, zero-based
      @param y: The y (vertical) coordinate, starts at the top, zero-based
      """
      idx = y * self.xSize + x
      ret = None
      
      pos = 0
      for cur, num in self.cmpData:
         pos = pos + num * len(cur)
         if pos > idx:
            tmp = list(cur)
            biggerList = num * tmp
            retIdx = idx - pos
            ret = biggerList[retIdx]
            break
            #ret = (num * list(cur))[idx - pos]
      return ret
   
   # ...........................
   def read(self, fn):
      """
      @summary: Reads in the compressed grid from a file
      @param fn: The filename to read from
      """
      with open(fn, 'rb') as f:
         method = struct.unpack('<B', f.read(1))[0]
         version = struct.unpack('<f', f.read(4))[0]
         self.xSize = struct.unpack('<l', f.read(4))[0]
         self.ySize = struct.unpack('<l', f.read(4))[0]

         tmpData = []
         
         # Read all data and compile list with keys
         tmp = f.read(1)
         while struct.unpack('<B', tmp)[0] != 0:
            cl = struct.unpack('<B', tmp)[0]
            num = struct.unpack('<B', f.read(1))[0]
            if num == 0:
               num = struct.unpack('<H', f.read(2))[0]
               if num == 0:
                  num = struct.unpack('<L', f.read(4))[0]
            tmpData.append((cl, num))
            tmp = f.read(1)
         
         # Read classes
         clDict = {}
         tmp = f.read(1)
         while tmp != "":
            clId = struct.unpack('<B', tmp)[0]
            nCats = struct.unpack('<B', f.read(1))[0]
            val = []
            for i in xrange(nCats):
               nVal = struct.unpack('<B', f.read(1))[0]
               val.append(nVal)
            clDict[clId] = val
            tmp = f.read(1)
         
         # Translate keys into classes
         self.cmpData = []
         
         for k, num in tmpData:
            self.cmpData.append((clDict[k], num))

   # ...........................
   def write(self, fn):
      """
      @summary: Writes out the compressed grid
      @param fn: The filename to write to
      """
      with open(fn, 'wb') as f:
         # Write method
         f.write(struct.pack('<B', self.METHOD))
         # Write version
         f.write(struct.pack('<f', self.VERSION))
         
         # Write x size and y size
         f.write(struct.pack('<l', self.xSize))
         f.write(struct.pack('<l', self.ySize))
         
         # Find all classes
         clDict = {}
         for cl, num in self.cmpData:
            if clDict.has_key(cl):
               clDict[cl] += 1
            else:
               clDict[cl] = 1
         
         sortedClasses = sorted(list(clDict.viewitems()), key=lambda k: k[1], reverse=True)
         
         #   Assign them ids
         clIds = {}
         i = 0
         for cl, num in sortedClasses:
            i += 1
            clIds[cl] = i

         #   Make sure that there are less than 256
         #      if not, break some of them down
         
         
         # Write data
         for cl, num in self.cmpData:
            # Write out class identifier
            f.write(struct.pack('<B', clIds[cl]))
            
            # Bit twidding based on run-length, if over threshold, we need to 
            #   write it out with more bytes.  We need to indicate this though
            #   so that it can be handled on read
            if num < 256:
               f.write(struct.pack('<B', num))
            elif num < 65536:
               # Write indicator to bump up size
               f.write(struct.pack('<B', 0))
               # Write number
               f.write(struct.pack('<H', num))
            else:
               # Write indicator to bump up twice
               f.write(struct.pack('<B', 0))
               f.write(struct.pack('<H', 0))
               # Write number
               f.write(struct.pack('<L', num))
         
         # Write separator
         f.write(struct.pack('<B', 0))
         
         # Write out all classes
         for k in clIds.keys():
            # Write id
            f.write(struct.pack('<B', clIds[k]))
            # Write number of items in key
            f.write(struct.pack('<B', len(list(k))))
            # Write items
            for i in list(k):
               f.write(struct.pack('<B', i))

# .............................................................................
if __name__ == "__main__":
   data = [
           [0, 0, 1, 0, 1, 1],
           [1, 2, 3, 2, 1, 0],
           [0, 3, 3, 2, 0, 0],
           [1, 3, 3, 1, 1, 0],
           [1, 3, 1, 1, 1, 1]
          ]
   
   mtx = Grid(griddedData=data)
   cmp = NormalRLECompressedGrid(grid=mtx)

   cmp.write('test2.bin')
   
   cmp2 = NormalRLECompressedGrid()
   cmp2.read('test2.bin')
   
   
   print cmp.query(1, 2)
   
   print cmp.decompress().data