"""
@summary: Module containing classes for Quadtree compression
"""
import struct

from matrix.matrix import Matrix

# .............................................................................
class QuadtreeHilbertCompressedMatrix(Matrix):
   # ...........................
   def __init__(self, mtx=None, rleThreshold=8):
      self.rleThreshold = rleThreshold
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
      
      self.data = quadtreeHilbertCompress(sqMtx, 2**self.rleThreshold)
   
   # ...........................
   def decompress(self):
      data = quadtreeHilbertDecompress(self.data, self.xSize)
      mtx = Matrix(data=data, xSize=self.xSize, ySize=self.ySize)
      return mtx
   
   # ...........................
   def getValue(self, x, y):
      # TODO
      minX = minY = 0
      maxX = maxY = self.xSize
      
      val = self.data
      
      while not isinstance(val, (int, list)):
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
      
      if isinstance(val, list):
         targetIdx = pointToHilbert(x, y, self.rleThreshold)
         idx = 0
         j = 0
         
         i, v = val[j]
         while j < len(val) and idx < targetIdx:
            idx = idx + i
            j += 1
            i, v = val[j]
         val = v
         
      return val

      
      return val
   
   # ...........................
   def readFile(self, filename):
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
      f.write(struct.pack('<B', 1))
      f.write(struct.pack('<B', data))
   elif isinstance(data, list):
      f.write(struct.pack('<B', 2))
      
      maxRL = 0
      count = 0
      for rl, _ in data:
         count += rl
         if rl > maxRL:
            maxRL = rl
      
      if maxRL < 256:
         mode = '<B'
         l = 1
      elif maxRL < 65536:
         mode = '<H'
         l = 2
      else:
         mode = '<L'
         l = 4
      
      order = 0
      while 2**order < count:
         order += 1
            
      # Write byte with starting bit and order
      firstByte = order + 128 * data[0][1]
      f.write(struct.pack('<B', firstByte))
      
      # Write how many bytes to read per entry
      f.write(struct.pack('<B', l))
      
      # Write entries
      for num, _ in data:
         f.write(struct.pack(mode, num))
   else:
      f.write(struct.pack('<B', 3))
      packData(data[1], f)
      packData(data[2], f)
      packData(data[3], f)
      packData(data[4], f)

# .............................................................................
def unpackData(f):
   itemType = struct.unpack('<B', f.read(1))[0]
   if itemType == 1:
      return struct.unpack('<B', f.read(1))[0]
   elif itemType == 2:
      firstByte = struct.unpack('<B', f.read(1))[0]
      if firstByte >= 128:
         bit = 1
         order = firstByte - 128
      else:
         bit = 0
         order = firstByte
      
      l = struct.unpack('<B', f.read(1))[0]
      if l == 1:
         mode = '<B'
      elif l == 2:
         mode = '<H'
      else:
         mode = '<L'
      
      ret = []
      count = 0
      while count < 2**order:
         num = struct.unpack(mode, f.read(l))[0]
         ret.append((num, bit))
         bit = -1 * bit + 1
         count += num
      return ret
   else:
      return {
              1: unpackData(f),
              2: unpackData(f),
              3: unpackData(f),
              4: unpackData(f)
             }

# .............................................................................
def quadtreeHilbertCompress(mtx, sizeThreshold):
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
      if numVals <= sizeThreshold:
         return hilbertCompress(mtx)
      else:
         h = len(mtx)
         return {
              1 : quadtreeHilbertCompress([mtx[i][:h/2] for i in range(h/2)], sizeThreshold),
              2 : quadtreeHilbertCompress([mtx[i][h/2:] for i in range(h/2)], sizeThreshold),
              3 : quadtreeHilbertCompress([mtx[i][:h/2] for i in range(h/2, h)], sizeThreshold),
              4 : quadtreeHilbertCompress([mtx[i][h/2:] for i in range(h/2, h)], sizeThreshold)
                }
# .............................................................................
def hilbertDecompress(data, order):
   lArray = []
   for count, val in data:
      lArray.extend(count*[val])

   mtx = []
   for y in xrange(2**order):
      row = []
      for x in xrange(2**order):
         row.append(lArray[pointToHilbert(x, y, order)])
      mtx.append(row)
         
   return mtx

# .............................................................................
def hilbertCompress(mtx):
   data = []
   xSize = len(mtx[0])
   ySize = len(mtx)

   order = 0
   while 2**order < xSize:
      order += 1

   lArray = 4**order * [0]
   
   for x in xrange(xSize):
      for y in xrange(ySize):
         idx = pointToHilbert(x, y, order)
         lArray[idx] = mtx[y][x]
   
   # Compress
   count = 1
   cur = lArray[0]
   for i in lArray[1:]:
      if i == cur:
         count += 1
      else:
         data.append((count, cur))
         cur = i
         count = 1
   data.append((count, cur))
   return data

# .............................................................................
def quadtreeHilbertDecompress(cmpMtx, sideLength):
   """
   @summary: Decompresses a quadtree compressed matrix
   @param cmpMtx: The compressed section of the matrix
   @param sideLength: The length of each side of this section
   """
   if isinstance(cmpMtx, int):
      return [[cmpMtx for i in xrange(sideLength)] for j in xrange(sideLength)]
   elif isinstance(cmpMtx, list):
      order = 0
      while 2**order < sideLength:
         order += 1
      return hilbertDecompress(cmpMtx, order)
   else:
      ret = []
      l = sideLength / 2
      topLeft = quadtreeHilbertDecompress(cmpMtx[1], l)
      topRight = quadtreeHilbertDecompress(cmpMtx[2], l)
      bottomLeft = quadtreeHilbertDecompress(cmpMtx[3], l)
      bottomRight = quadtreeHilbertDecompress(cmpMtx[4], l)
      for i in xrange(l):
         ret.append(topLeft[i] + topRight[i])
      for i in xrange(l):
         ret.append(bottomLeft[i] + bottomRight[i])
      return ret

# .............................................................................
hilbertMap = { 
              'a': {
                    (0, 0): (0, 'd'), 
                    (0, 1): (1, 'a'), 
                    (1, 0): (3, 'b'), 
                    (1, 1): (2, 'a')
                   }, 
              'b': {
                    (0, 0): (2, 'b'), 
                    (0, 1): (1, 'b'), 
                    (1, 0): (3, 'a'), 
                    (1, 1): (0, 'c')
                   }, 
              'c': {
                    (0, 0): (2, 'c'), 
                    (0, 1): (3, 'd'), 
                    (1, 0): (1, 'c'), 
                    (1, 1): (0, 'b')
                   }, 
              'd': {
                    (0, 0): (0, 'a'), 
                    (0, 1): (3, 'c'), 
                    (1, 0): (1, 'd'), 
                    (1, 1): (2, 'd')
                   } 
              }

# .............................................................................
def pointToHilbert(x, y, order): 
   current_square = 'a' 
   position = 0 
   for i in range(order - 1, -1, -1): 
      position <<= 2 
      quad_x = 1 if x & (1 << i) else 0 
      quad_y = 1 if y & (1 << i) else 0 
      quad_position, current_square = hilbertMap[current_square][(quad_x, quad_y)] 
      position |= quad_position 
   return position

   
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
   cmpMtx = QuadtreeHilbertCompressedMatrix(rleThreshold=2)
   cmpMtx.compress(mtx)
   cmpMtx.writeFile(fn)
   #print cmpMtx.data
   #mtx2 = cmpMtx.decompress()
   #for row in mtx2.data:
   #   print row
   mtx2 = QuadtreeHilbertCompressedMatrix()
   mtx2.readFile(fn)
   print mtx2.data
#   mtx2.readFile(fn)
   mtx3 = mtx2.decompress()
   for row in mtx3.data:
      print row
   
   m = []
   for y in xrange(cmpMtx.ySize):
      r = []
      for x in xrange(cmpMtx.xSize):
         r.append(cmpMtx.getValue(x, y))
      m.append(r)

   for r in m:
      print r
      
if __name__ == "__main__":
   test()
