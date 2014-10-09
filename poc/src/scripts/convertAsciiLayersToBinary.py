"""
@summary: Script to convert all of the layer ASCII data to binary data files
@author: CJ Grady
"""
import os
from tools.ascii import convertAscToMatrix

inDir = '/data/geo716/final project/data/ascs'
outDir = '/data/geo716/final project/data/bins'

for i in xrange(1000):
   #i = 1
   print i
   ascFn = os.path.join(inDir, '%s.asc' % i)
   binFn = os.path.join(outDir, '%s.bin' % i)
   mtx = convertAscToMatrix(ascFn, threshold=10)
   mtx.writeFile(binFn)
  
