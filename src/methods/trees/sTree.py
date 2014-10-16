"""
@summary: Module containing a class to perform S-Tree compression
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
class STreeCompressedGrid(_CompressedGrid):
   """
   @summary: This class compresses a grid using an S-Tree.  This basically 
                divides a grid in half at each step.
   """
   
   METHOD = 4
   VERSION = 2.0
   
   # ...........................
   def __init__(self, mtx=None):
      if mtx is not None:
         self.compress(mtx)
      else:
         self.data = []
         self.xSize = None
         self.ySize = None
         self.order = 0
         
   # ...........................
   def compress(self, mtx):
      self.data = []
      self.xSize = mtx.xSize
      self.ySize = mtx.ySize
      
      assert len(mtx.data) == mtx.ySize
      assert len(mtx.data[0]) == mtx.xSize
      
      # Determine order
      self.xOrder = self.yOrder = 0
      
      while 2**self.xOrder < self.xSize:
         self.xOrder += 1
         
      while 2**self.yOrder < self.ySize:
         self.yOrder += 1
      
      self.order = 0
      while 2**self.order < max([self.xSize, self.ySize]):
         self.order += 1
      
      paddedMtx = []
      for row in mtx.data:
         paddedMtx.append(row + (2**self.xOrder - len(row)) * [0])
      for i in xrange(len(paddedMtx), 2**self.yOrder):
         paddedMtx.append(2**self.xOrder * [0])
      
      self.data = streeCompress(paddedMtx, self.xOrder, self.yOrder)
   
   # ...........................
   def decompress(self):
      data = streeDecompress(self.data, self.xSize, self.ySize)
      mtx = Grid(griddedData=data)
      return mtx
   
   # ...........................
   def query(self, x, y):
      minX = minY = 0
      maxX = self.xSize
      maxY = self.ySize
      
      val = self.data
      
      while not isinstance(val, int):
         
         if val['splitX']:
            cmpX = (maxX+minX) / 2
            if x < cmpX:
               maxX = cmpX
               val = val[1]
            else:
               minX = cmpX
               val = val[2]
         else:
            cmpY = (maxY+minY) / 2
            if y < cmpY:
               maxY = cmpY
               val = val[1]
            else:
               minY = cmpY
               val = val[2]
   
      return val
   
   # ...........................
   def read(self, filename):
      with open(filename, 'rb') as f:
         self.xSize = struct.unpack('<l', f.read(4))[0]
         self.ySize = struct.unpack('<l', f.read(4))[0]

         self.data = unpackData(f)
   
   # ...........................
   def write(self, filename):
      with open(filename, 'wb') as f:
         # write version
         f.write(struct.pack('<l', self.xSize))
         f.write(struct.pack('<l', self.ySize))
         packData(self.data, f)

# .............................................................................
def packData(data, f):
   if isinstance(data, int):
      f.write(struct.pack('<?', True))
      f.write(struct.pack('<B', data))
   else:
      f.write(struct.pack('<?', False))
      f.write(struct.pack('<?', data['splitX']))
      packData(data[1], f)
      packData(data[2], f)

# .............................................................................
def unpackData(f):
   isInt = struct.unpack('<?', f.read(1))[0]
   if isInt:
      return struct.unpack('<B', f.read(1))[0]
   else:
      return {
              'splitX' : struct.unpack('<?', f.read(1))[0],
              1: unpackData(f),
              2: unpackData(f)
             }

# .............................................................................
def streeCompress(mtx, xOrder, yOrder):
   """
   @summary: Performs s-tree compression (approximate)
   @param mtx: List of lists
   @param xOrder: 2**xOrder elements in each row
   @param yOrder: 2**yOrder rows
   """
   minV = min([min(r) for r in mtx])
   maxV = max([max(r) for r in mtx])

   if minV == maxV:
      return minV
   else:
      if xOrder > yOrder:
         return {
                 'splitX' : True,
                 1: streeCompress([r[:2**(xOrder-1)] for r in mtx], xOrder-1, yOrder),
                 2: streeCompress([r[2**(xOrder-1):] for r in mtx], xOrder-1, yOrder)
                }
      else:
         return {
                 'splitX' : False,
                 1: streeCompress(mtx[:2**(yOrder-1)], xOrder, yOrder-1),
                 2: streeCompress(mtx[2**(yOrder-1):], xOrder, yOrder-1)
                }

# .............................................................................
def streeDecompress(cmpMtx, xSize, ySize):
   """
   @summary: Decompresses an s-tree compressed matrix
   @param cmpMtx: The compressed section of the matrix
   """
   if isinstance(cmpMtx, int):
      return [[cmpMtx for i in xrange(xSize)] for j in xrange(ySize)]
   else:
      ret = []
      if cmpMtx['splitX']:
         leftSide = streeDecompress(cmpMtx[1], xSize/2, ySize)
         rightSide = streeDecompress(cmpMtx[2], xSize/2, ySize)
         for i in xrange(ySize):
            ret.append(leftSide[i] + rightSide[i])
      else:
         ret = streeDecompress(cmpMtx[1], xSize, ySize/2)
         ret.extend(streeDecompress(cmpMtx[2], xSize, ySize/2))
      return ret
   
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
          [1, 0, 0, 1, 0, 1, 1, 1, 4, 4, 4, 4, 4, 4, 4, 4], 
          [0, 1, 0, 1, 0, 1, 1, 1, 4, 4, 4, 4, 4, 4, 4, 4], 
          [0, 1, 1, 1, 1, 1, 1, 0, 4, 4, 4, 4, 4, 4, 4, 4], 
          [1, 0, 1, 1, 0, 0, 1, 1, 4, 4, 4, 4, 4, 4, 4, 4], 
          [1, 1, 0, 1, 1, 1, 1, 1, 4, 4, 4, 4, 4, 4, 4, 4], 
          [0, 0, 0, 1, 1, 1, 1, 0, 4, 4, 4, 4, 4, 4, 4, 4], 
          [0, 0, 0, 0, 0, 1, 0, 1, 4, 4, 4, 4, 4, 4, 4, 4], 
          [1, 0, 0, 0, 1, 0, 1, 1, 4, 4, 4, 4, 4, 4, 4, 4]
         ]
   mtx = Grid(data)
   cmpMtx = STreeCompressedGrid()
   cmpMtx.compress(mtx)
   
   g2 = cmpMtx.decompress()
   #cmpMtx.writeFile(fn)
   print cmpMtx.data
   print g2.data
#    fn = '/home/cjgrady/test.qtree'
#    mtx = Matrix(data=data)
#    cmpMtx = STreeCompressedMatrix()
#    cmpMtx.compress(mtx)
#    cmpMtx.writeFile(fn)
#    #print cmpMtx.data
#    #mtx2 = cmpMtx.decompress()
#    mtx2 = STreeCompressedMatrix()
#    mtx2.readFile(fn)
#    mtx3 = mtx2.decompress()
#    for row in mtx3.data:
#       print row
#    print cmpMtx.getValue(1, 0)
#      
      
if __name__ == "__main__":
   test()
