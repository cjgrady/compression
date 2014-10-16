"""
@summary: This module contains a Quadtree variant class that uses a Quadtree
             compression technique until a threshold is met, then switches to
             a run-length encoding technique to fill in details
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
from methods.combo.comboBase import _ComboCompressionBaseClass

from methods.rle.hilbert import HilbertRLECompressedGrid

# .............................................................................
class QuadtreeComboCompressedGrid(_ComboCompressionBaseClass):
   """
   @summary: This class provides a two-stage approach for compressing a grid by
                first using a Quadtree and then switching to a run-length
                encoding mechanism once a threshold is reached
   """
   METHOD = 5
   VERSION = 2.0
   
   def _initialize(self):
      self.cmpData = {}
      self.xSize = None
      self.ySize = None
      self.order = 0

   # ...........................
   def compress(self, mtx):
      self.cmpData = {}
      self.xSize = mtx.xSize
      self.ySize = mtx.ySize
      
      assert len(mtx.data) == mtx.ySize
      assert len(mtx.data[0]) == mtx.xSize
      
      # Determine order
      self.order = 0
      while 2**self.order < max([self.xSize, self.ySize]):
         self.order += 1
      
      sqMtx = []
      for row in mtx.data:
         sqMtx.append(row + (2**self.order - len(row)) * [0])
      for i in xrange(len(sqMtx), 2**self.order):
         sqMtx.append(2**self.order * [0])
      
      self.cmpData = quadtreeComboCompress(sqMtx, self.rleMethod, self.threshold)
   
   # ...........................
   def decompress(self):
      data = quadtreeComboDecompress(self.cmpData, self.xSize)
      mtx = Grid(griddedData=data)
      return mtx

   # ...........................
   def write(self, fn):
      """
      @summary: Writes out the compressed grid
      @param fn: The filename to write to
      """
       
      # ...........................
      # Find categories
      
      # Find all classes
      clDict = {}
       
      def findCategories(d):
         cats = {}
         if isinstance(d, int):
            return {}
         elif isinstance(d, _CompressedGrid):
            # Do stuff here
            for cl, _ in d.cmpData:
               if cats.has_key(cl):
                  cats[cl] += 1
               else:
                  cats[cl] = 1
         else:
            dicts = [findCategories(d[1]),
                     findCategories(d[2]),
                     findCategories(d[3]),
                     findCategories(d[4])]
            for catDict in dicts:
               for k in catDict.keys():
                  if cats.has_key(k):
                     cats[k] += catDict[k]
                  else:
                     cats[k] = catDict[k]
         return cats
          
      clDict = findCategories(self.cmpData)
       
      sortedClasses = sorted(list(clDict.viewitems()), key=lambda k: k[1], reverse=True)
       
      #   Assign them ids
      clIds = {}
      i = 0
      for cl, _ in sortedClasses:
         i += 1
         clIds[cl] = i
         
      # Write out results
      with open(fn, 'wb') as f:
         # Write method
         f.write(struct.pack('<B', self.METHOD))
         # Write version
         f.write(struct.pack('<f', self.VERSION))
         
         # RLE method and version
         f.write(struct.pack('<B', self.rleMethod.METHOD))
         f.write(struct.pack('<f', self.rleMethod.VERSION))
         
         # Write x size and y size
         f.write(struct.pack('<l', self.xSize))
         f.write(struct.pack('<l', self.ySize))
         
         
         def packData(data, f):
            if isinstance(data, int):
               f.write(struct.pack('<B', 1))
               f.write(struct.pack('<B', data))
            elif isinstance(data, _CompressedGrid):
               f.write(struct.pack('<B', 2))

               # Might need to write x and y sizes
               # Write data
               for cl, num in data.cmpData:
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
      
            else:
               f.write(struct.pack('<B', 3))
               packData(data[1], f)
               packData(data[2], f)
               packData(data[3], f)
               packData(data[4], f)
      
         # Pack data
         packData(self.cmpData, f)
      
         # Write out categories
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
def quadtreeComboCompress(mtx, method, sizeThreshold):
   """
   @summary: Performs quadtree compression
   @param mtx: List of lists, assumed to be square
   """
   minV = min([min(r) for r in mtx])
   maxV = max([max(r) for r in mtx])

   if minV == maxV:
      return minV
   else:
      numVals = len(mtx)**2
      if numVals <= sizeThreshold:
         return method(grid=Grid(griddedData=mtx))
      else:
         h = len(mtx)
         return {
              1 : quadtreeComboCompress([mtx[i][:h/2] for i in range(h/2)], method, sizeThreshold),
              2 : quadtreeComboCompress([mtx[i][h/2:] for i in range(h/2)], method, sizeThreshold),
              3 : quadtreeComboCompress([mtx[i][:h/2] for i in range(h/2, h)], method, sizeThreshold),
              4 : quadtreeComboCompress([mtx[i][h/2:] for i in range(h/2, h)], method, sizeThreshold)
                }
   
# .............................................................................
def quadtreeComboDecompress(cmpMtx, sideLength):
   """
   @summary: Decompresses a quadtree compressed matrix
   @param cmpMtx: The compressed section of the matrix
   @param sideLength: The length of each side of this section
   """
   if isinstance(cmpMtx, int):
      return [[cmpMtx for i in xrange(sideLength)] for j in xrange(sideLength)]
   elif isinstance(cmpMtx, _CompressedGrid):
      return cmpMtx.decompress().data
   else:
      ret = []
      l = sideLength / 2
      topLeft = quadtreeComboDecompress(cmpMtx[1], l)
      topRight = quadtreeComboDecompress(cmpMtx[2], l)
      bottomLeft = quadtreeComboDecompress(cmpMtx[3], l)
      bottomRight = quadtreeComboDecompress(cmpMtx[4], l)
      for i in xrange(l):
         ret.append(topLeft[i] + topRight[i])
      for i in xrange(l):
         ret.append(bottomLeft[i] + bottomRight[i])
      return ret


   
   
   
   
   
#    # ...........................
#    def query(self, x, y):
#       minX = minY = 0
#       maxX = maxY = self.xSize
#       
#       val = self.data
#       
#       while not isinstance(val, int):
#          cmpX = (maxX+minX) / 2
#          cmpY = (maxY+minY) / 2
#          if x < cmpX and y < cmpY:
#             maxX = cmpX
#             maxY = cmpY
#             val = val[1]
#          elif x < cmpX:
#             maxX = cmpX
#             minY = cmpY
#             val = val[3]
#          elif y < cmpY:
#             minX = cmpX
#             maxY = cmpY
#             val = val[2]
#          else:
#             minX = cmpX
#             minY = cmpY
#             val = val[4]
#    
#       return val
#    
#    # ...........................
#    def read(self, filename):
#       # TODO
#       with open(filename, 'rb') as f:
#          self.xSize = struct.unpack('<l', f.read(4))[0]
#          self.ySize = struct.unpack('<l', f.read(4))[0]
# 
#          self.data = unpackData(f)
#    
#    # ...........................
#    def write(self, filename):
#       with open(filename, 'wb') as f:
#          # write version
#          f.write(struct.pack('<l', self.xSize))
#          f.write(struct.pack('<l', self.ySize))
#          packData(self.data, f)

# .............................................................................
def packData(data, f):
   if isinstance(data, int):
      f.write(struct.pack('<?', True))
      f.write(struct.pack('<B', data))
   else:
      f.write(struct.pack('<?', False))
      packData(data[1], f)
      packData(data[2], f)
      packData(data[3], f)
      packData(data[4], f)

# .............................................................................
def unpackData(f):
   isInt = struct.unpack('<?', f.read(1))[0]
   if isInt:
      return struct.unpack('<B', f.read(1))[0]
   else:
      return {
              1: unpackData(f),
              2: unpackData(f),
              3: unpackData(f),
              4: unpackData(f)
             }

   
# .............................................................................
def test():
   data = [
          [0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0], 
          [1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0], 
          [1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0], 
          [0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1], 
          [0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1], 
          [0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0], 
          [1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1], 
          [0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1], 
          [1, 0, 0, 1, 0, 1, 1, 1, 3, 0, 1, 0, 0, 0, 1, 1], 
          [0, 1, 0, 1, 0, 1, 1, 1, 4, 1, 0, 0, 0, 1, 0, 1], 
          [0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1], 
          [1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1], 
          [1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1], 
          [0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1], 
          [2, 2, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1], 
          [2, 2, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1]
         ]
   mtx = Grid(data)
   cmpMtx = QuadtreeComboCompressedGrid(HilbertRLECompressedGrid, orderThreshold=2)
   cmpMtx.compress(mtx)
   #cmpMtx.writeFile(fn)
   print cmpMtx.cmpData
   
   g2 = cmpMtx.decompress()
   print
   for row in g2.data:
      print row
   #print g2.data
#   mtx2 = QuadtreeCompressedMatrix()
#   mtx2.readFile(fn)
#   mtx3 = mtx2.decompress()
#   for row in mtx3.data:
#      print row
#      
      
if __name__ == "__main__":
   test()

