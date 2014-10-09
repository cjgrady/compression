"""
@summary: Module containing Matrix base class
@author: CJ Grady
"""
import struct

# .............................................................................
class Matrix(object):
   # ...........................
   def __init__(self, data=None, xSize=None, ySize=None):
      if data is not None:
         self.data = data
         self.xSize = xSize if xSize is not None else len(self.data[0])
         self.ySize = ySize if ySize is not None else len(self.data)
         assert self.xSize == len(self.data[0])
         assert self.ySize == len(self.data)
      else:
         self.data = []
         self.xSize = None
         self.ySize = None

   # ...........................
   def getValue(self, x, y):
      return self.data[y][x]
   
   # ...........................
   def readFile(self, filename):
      with open(filename, 'rb') as f:
         self.xSize = struct.unpack('<l', f.read(4))[0]
         self.ySize = struct.unpack('<l', f.read(4))[0]
         self.data = []
         
         i = 7
         x = 0
         y = 0
         row = []
         # Read data
         byte = f.read(1)
         while byte != "":
            num = struct.unpack('<B', byte)[0]
            for i in range(7, -1, -1):
               val = num / 2**i
               if val:
                  num = num - 2**i
               row.append(val)
               x += 1
               if x == self.xSize:
                  y += 1
                  x = 0
                  self.data.append(row)
                  row = []
                  
            byte = f.read(1)
   
   # ...........................
   def writeFile(self, filename):
      with open(filename, 'wb') as f:
         # write version
         f.write(struct.pack('<l', self.xSize))
         f.write(struct.pack('<l', self.ySize))
         # write data
         lArray = []
         for row in self.data:
            lArray.extend(row)
         i = 7
         num = 0
         for j in lArray:
            num = num + j * 2**i
            if i == 0:
               f.write(struct.pack('<B', num))
               num = 0
               i = 7
            else:
               i -= 1
         if i < 7:
            f.write(struct.pack('<B', num))

if __name__ == "__main__":
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
   fn = '/home/cjgrady/cjTest.bin'
   
   mtx = Matrix(data=data)
   mtx.writeFile(fn)

   mtx2 = Matrix()
   mtx2.readFile(fn)
   
   for row in mtx2.data:
      print row
      #print mtx2.data
   
