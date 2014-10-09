"""
@summary: Module containing classes for normal run length encoding
"""
import struct

from matrix.matrix import Matrix

# .............................................................................
class NormalRLECompressedMatrix(Matrix):
   # ...........................
   def __init__(self, mtx=None):
      if mtx is not None:
         self.compress(mtx)
      else:
         self.data = []
         self.xSize = None
         self.ySize = None
         
   # ...........................
   def compress(self, mtx):
      self.data = []
      self.xSize = mtx.xSize
      self.ySize = mtx.ySize
      
      lArray = []
      for row in mtx.data:
         lArray.extend(row)
      
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
      data = [lArray[i:i+self.xSize] for i in xrange(0, self.xSize*self.ySize, self.xSize)]
      mtx = Matrix(data=data, xSize=self.xSize, ySize=self.ySize)
      return mtx
   
   # ...........................
   def getValue(self, x, y):
      targetIdx = y*self.xSize + x
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
   cmpMtx = NormalRLECompressedMatrix()
   cmpMtx.compress(mtx)
   mtx2 = cmpMtx.decompress()

   cmpMtx.writeFile('/home/cjgrady/cjTestRLE.bin')
   mtx3 = NormalRLECompressedMatrix()
   mtx3.readFile('/home/cjgrady/cjTestRLE.bin')
   
      
      
if __name__ == "__main__":
   test()
   