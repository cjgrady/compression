"""
@summary: Module containing classes for Quadtree compression
"""
import struct

from matrix.matrix import Matrix

# .............................................................................
class QuadtreeCompressedMatrix(Matrix):
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
      self.order = 0
      while 2**self.order < max([self.xSize, self.ySize]):
         self.order += 1
      
      sqMtx = []
      for row in mtx.data:
         sqMtx.append(row + (2**self.order - len(row)) * [0])
      for i in xrange(len(sqMtx), 2**self.order):
         sqMtx.append(2**self.order * [0])
      
      self.data = quadtreeCompress(sqMtx)
   
   # ...........................
   def decompress(self):
      data = quadtreeDecompress(self.data, self.xSize)
      mtx = Matrix(data=data, xSize=self.xSize, ySize=self.ySize)
      return mtx
   
   # ...........................
   def getValue(self, x, y):
      minX = minY = 0
      maxX = maxY = self.xSize
      
      val = self.data
      
      while not isinstance(val, int):
         cmpX = (maxX+minX) / 2
         cmpY = (maxY+minY) / 2
         if x < cmpX and y < cmpY:
            maxX = cmpX
            maxY = cmpY
            val = val[1]
         elif x < cmpX:
            maxX = cmpX
            minY = cmpY
            val = val[3]
         elif y < cmpY:
            minX = cmpX
            maxY = cmpY
            val = val[2]
         else:
            minX = cmpX
            minY = cmpY
            val = val[4]
   
      return val
   
   # ...........................
   def readFile(self, filename):
      # TODO
      with open(filename, 'rb') as f:
         self.xSize = struct.unpack('<l', f.read(4))[0]
         self.ySize = struct.unpack('<l', f.read(4))[0]

         self.data = unpackData(f)
   
   # ...........................
   def writeFile(self, filename):
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
def quadtreeCompress(mtx):
   """
   @summary: Performs quadtree compression
   @param mtx: List of lists, assumed to be square
   """
   sumOfMtx = sum([sum(r) for r in mtx])
   numVals = len(mtx)**2
   if sumOfMtx == 0:  # All zeros
      return 0
   elif sumOfMtx == numVals: # All ones
      return 1
   else:
      h = len(mtx)
      return {
              1 : quadtreeCompress([mtx[i][:h/2] for i in range(h/2)]),
              2 : quadtreeCompress([mtx[i][h/2:] for i in range(h/2)]),
              3 : quadtreeCompress([mtx[i][:h/2] for i in range(h/2, h)]),
              4 : quadtreeCompress([mtx[i][h/2:] for i in range(h/2, h)])
             }

# .............................................................................
def quadtreeDecompress(cmpMtx, sideLength):
   """
   @summary: Decompresses a quadtree compressed matrix
   @param cmpMtx: The compressed section of the matrix
   @param sideLength: The length of each side of this section
   """
   if isinstance(cmpMtx, int):
      return [[cmpMtx for i in xrange(sideLength)] for j in xrange(sideLength)]
   else:
      ret = []
      l = sideLength / 2
      topLeft = quadtreeDecompress(cmpMtx[1], l)
      topRight = quadtreeDecompress(cmpMtx[2], l)
      bottomLeft = quadtreeDecompress(cmpMtx[3], l)
      bottomRight = quadtreeDecompress(cmpMtx[4], l)
      for i in xrange(l):
         ret.append(topLeft[i] + topRight[i])
      for i in xrange(l):
         ret.append(bottomLeft[i] + bottomRight[i])
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
          [1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1], 
          [0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1], 
          [0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1], 
          [1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1], 
          [1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1], 
          [0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1], 
          [0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1], 
          [1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1]
         ]
   fn = '/home/cjgrady/test.qtree'
   mtx = Matrix(data=data)
   cmpMtx = QuadtreeCompressedMatrix()
   cmpMtx.compress(mtx)
   cmpMtx.writeFile(fn)
   print cmpMtx.data
#   mtx2 = QuadtreeCompressedMatrix()
#   mtx2.readFile(fn)
#   mtx3 = mtx2.decompress()
#   for row in mtx3.data:
#      print row
#      
      
if __name__ == "__main__":
   test()
