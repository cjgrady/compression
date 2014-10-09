"""
@summary: Module containing classes for run length encoding using a Hilbert
             curve
"""
import struct

from matrix.matrix import Matrix

# .............................................................................
class HilbertRLECompressedMatrix(Matrix):
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
      
      # Initialize lArray
      lArray = 4**self.order * [0]
      
      # go through matrix and add items to array
      for x in xrange(self.xSize):
         for y in xrange(self.ySize):
            idx = pointToHilbert(x, y, self.order)
            try:
               lArray[idx] = mtx.data[y][x]
            except:
               print "Max x:", self.xSize
               print "Max y:", self.ySize
               print "idx:", idx
               print "x:", x
               print "y:", y
               print "len larray:", len(lArray)
               print "len data:", len(mtx.data)
               print "len data[y]:", len(mtx.data[y])
      
      # Compress
      count = 1
      cur = lArray[0]
      for i in lArray[1:]:
         if i == cur:
            count += 1
         else:
            self.data.append((count, cur))
            cur = i
            count = 1
      self.data.append((count, cur))
   
   # ...........................
   def decompress(self):
      lArray = []
      for count, val in self.data:
         lArray.extend(count*[val])
         
      data = []
      for y in xrange(self.ySize):
         row = []
         for x in xrange(self.xSize):
            row.append(lArray[pointToHilbert(x, y, self.order)])
         data.append(row)
      
      mtx = Matrix(data=data, xSize=self.xSize, ySize=self.ySize)
      return mtx
   
   # ...........................
   def getValue(self, x, y):
      targetIdx = pointToHilbert(x, y, self.order)
      idx = 0
      j = 0
      
      i, val = self.data[j]
      while j < len(self.data) and idx < targetIdx:
         idx = idx + i
         j += 1
         i, val = self.data[j]
      
      return val
   
   # ...........................
   def readFile(self, filename):
      with open(filename, 'rb') as f:
         self.xSize = struct.unpack('<l', f.read(4))[0]
         self.ySize = struct.unpack('<l', f.read(4))[0]
         bit = struct.unpack('<B', f.read(1))[0]
         l = struct.unpack('<B', f.read(1))[0]
         self.data = []

         if l == 1:
            mode = '<B'
         elif l == 2:
            mode = '<H'
         else:
            mode = '<L'
         
         # Read data
         byte = f.read(l)
         while byte != "":
            self.data.append((struct.unpack(mode, byte)[0], bit))
            bit = (bit -1)**2
            byte = f.read(l)
   
   # ...........................
   def writeFile(self, filename):
      with open(filename, 'wb') as f:
         # write version
         f.write(struct.pack('<l', self.xSize))
         f.write(struct.pack('<l', self.ySize))
         f.write(struct.pack('<B', self.data[0][1]))

         # Find largest run
         maxRL = 0
         for rl, _ in self.data:
            if rl > maxRL:
               maxRL = rl
         
         # Determine mode to use
         if maxRL < 256:
            mode = '<B'
            l = 1
         elif maxRL < 65536:
            mode = '<H'
            l = 2
         else:
            mode = '<L'
            l = 4
            
         # Write how many bits to read per item
         f.write(struct.pack('<B', l))
         
         # write data
         for count, _ in self.data:
            f.write(struct.pack(mode, count))

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
   mtx = Matrix(data=data)
   cmpMtx = HilbertRLECompressedMatrix()
   cmpMtx.compress(mtx)
   #print cmpMtx.data
   mtx2 = cmpMtx.decompress()
   for row in mtx2.data:
      print row
      
      
if __name__ == "__main__":
   test()
