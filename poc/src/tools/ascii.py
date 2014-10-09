"""
@summary: Module containing ascii tools
"""

from matrix.matrix import Matrix
from rle.normal import NormalRLECompressedMatrix

def convertAscToMatrix(asciiFn, threshold=1):
   # Assume ncols, nrows, then 3 other metadata lines before data
   f = open(asciiFn, 'r')
   lines = f.readlines()
   f.close()
   xSize = int(lines[0].split('ncols')[1].strip())
   ySize = int(lines[1].split('nrows')[1].strip())
   # 2 - xll
   # 3 - yll
   # 4 - cell size
   # 5 - no data
   NODATA_VALUE = -1
   #NODATA_VALUE = int(lines[5].split('NODATA_value')[1].strip())
   data = []
   for line in lines[5:]:
      row = []
      vals = line.split(' ')
      for j in vals:
         if len(j) > 0:
            i = float(j)
            if i < 1.0:
               i = i * 100
            i = int(i)
            
            #for i in [int(j) for j in line.split('  ')]:
            if i >= threshold and i != NODATA_VALUE:
               row.append(1)
            else:
               row.append(0)
      data.append(row)
   mtx = Matrix(data=data, xSize=xSize, ySize=ySize)
   return mtx

if __name__ == "__main__":
   fn = '/home/cjgrady/surfaces2/miPos0.9970.asc'
   mtx = convertAscToMatrix(fn, threshold=2)
   for row in mtx.data:
      print row
   mtx.writeFile('/home/cjgrady/test1.bin')
      
   mtx2 = NormalRLECompressedMatrix()
   mtx2.compress(mtx)
   
   print mtx2.data
   
   mtx2.writeFile('/home/cjgrady/test.rlebin')
   