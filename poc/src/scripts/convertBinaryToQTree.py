"""
@summary: Module containing script to convert binary rasters to compressed 
             rasters via quadtree encoding
"""
import os
from matrix.matrix import Matrix
from tree.quadTree import QuadtreeCompressedMatrix

inDir = '/data/geo716/final project/data/bins'
outDir = '/data/geo716/final project/data/qtrees'

for i in xrange(1000):
   #i = 0
   print i
   binFn = os.path.join(inDir, '%s.bin' % i)
   rleFn = os.path.join(outDir, '%s.qtree' % i)
   mtx = Matrix()
   mtx.readFile(binFn)
   mtx2 = QuadtreeCompressedMatrix()
   mtx2.compress(mtx)
   mtx2.writeFile(rleFn)
  