"""
@summary: Module containing script to convert binary rasters to compressed 
             rasters via combo quad tree hilbert curve encoding
"""
import os
from matrix.matrix import Matrix
from combo.quadTreeHilbert import QuadtreeHilbertCompressedMatrix

inDir = '/data/geo716/final project/data/bins'
outDir = '/data/geo716/final project/data/combos'

orders = [6, 8, 10]

for i in xrange(1000):
   binFn = os.path.join(inDir, '%s.bin' % i)
   mtx = Matrix()
   mtx.readFile(binFn)
   for order in orders:
      print i, '-', order
      cmpFn = os.path.join('%s%s' % (outDir, order), '%s.combo%s' % (i, order))
      
      mtx2 = QuadtreeHilbertCompressedMatrix(rleThreshold=order)
      mtx2.compress(mtx)
      mtx2.writeFile(cmpFn)
