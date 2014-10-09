"""
@summary: Module containing script to convert binary rasters to compressed 
             rasters via hilbert curve run length encoding
"""
import os
from matrix.matrix import Matrix
from rle.hilbert import HilbertRLECompressedMatrix

inDir = '/data/geo716/final project/data/bins'
outDir = '/data/geo716/final project/data/hilberts'

for i in xrange(1000):
   #i = 0
   print i
   binFn = os.path.join(inDir, '%s.bin' % i)
   rleFn = os.path.join(outDir, '%s.hilb' % i)
   mtx = Matrix()
   mtx.readFile(binFn)
   mtx2 = HilbertRLECompressedMatrix()
   mtx2.compress(mtx)
   mtx2.writeFile(rleFn)
  
