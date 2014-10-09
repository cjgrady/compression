"""
@summary: Module containing script to convert binary rasters to compressed 
             rasters via normal run length encoding
"""
import os
from matrix.matrix import Matrix
from rle.normal import NormalRLECompressedMatrix

inDir = '/data/geo716/final project/data/bins'
outDir = '/data/geo716/final project/data/rles'

for i in xrange(1000):
   #i = 0
   print i
   binFn = os.path.join(inDir, '%s.bin' % i)
   rleFn = os.path.join(outDir, '%s.rle' % i)
   mtx = Matrix()
   mtx.readFile(binFn)
   mtx2 = NormalRLECompressedMatrix()
   mtx2.compress(mtx)
   mtx2.writeFile(rleFn)
  
