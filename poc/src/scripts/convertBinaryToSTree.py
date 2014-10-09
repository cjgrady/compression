"""
@summary: Module containing script to convert binary rasters to compressed 
             rasters via stree encoding
"""
import os
from matrix.matrix import Matrix
from tree.sTree import STreeCompressedMatrix

inDir = '/data/geo716/final project/data/bins'
outDir = '/data/geo716/final project/data/strees'

for i in xrange(1000):
   #i = 0
   print i
   binFn = os.path.join(inDir, '%s.bin' % i)
   rleFn = os.path.join(outDir, '%s.stree' % i)
   mtx = Matrix()
   mtx.readFile(binFn)
   mtx2 = STreeCompressedMatrix()
   mtx2.compress(mtx)
   mtx2.writeFile(rleFn)
  