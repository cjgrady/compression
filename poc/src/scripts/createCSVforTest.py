"""
@summary: Creates a CSV file with all of the measurements and stats collected
@author: CJ Grady
"""
import csv
import glob
import os

from combo.quadTreeHilbert import QuadtreeHilbertCompressedMatrix
from rle.hilbert import HilbertRLECompressedMatrix
from rle.normal import NormalRLECompressedMatrix
from tree.quadTree import QuadtreeCompressedMatrix
from tree.sTree import STreeCompressedMatrix

from tools.ascii import convertAscToMatrix

# .............................................................................
def getFilename(baseDir, subDir, name, ext):
   return os.path.join(baseDir, subDir, '%s.%s' % (name, ext))

# .............................................................................
def getFileSize(fn):
   return os.path.getsize(fn)

# .............................................................................
def getNumPresencesAndVariance(mtx):
   numCells = len(mtx) * len(mtx[0])
   numPresences = sum([sum(r) for r in mtx])
   mean = 1.0 * numPresences / numCells
   a = 0.0
   for row in mtx:
      for col in row:
         a += (col - mean)**2
   variance = 1.0 * a / numCells
   return numPresences, variance
   

# Find all files to check

DATA_DIR = "/data/geo716/final project/data/tests/"
ASC = 'ascs'
TIFF = 'tiffs'
BIN = 'bins'
RLE = 'rles'
HILBERT = 'hilberts'
QTREE = 'qtrees'
STREE = 'strees'
COMBO6 = 'combo6'
COMBO8 = 'combo8'
COMBO10 = 'combo10'


OUT_CSV = "/data/geo716/final project/data/testsOutput.csv"

fns = glob.glob('%s/%s/*' % (DATA_DIR, ASC))

headerRow = [
             'Number of Presences', 'Morans I', 'Variance', 'ASCII Size', 
             'Tiff Size', 'Binary Size', 'RLE Size', 'Hilbert Size', 
             'Quadtree Size', 'S-Tree Size', 'Combo 6 Size', 'Combo 8 Size', 
             'Combo 10 Size'
            ]

with open(OUT_CSV, 'w') as csvOutFile:
   writer = csv.writer(csvOutFile)
   writer.writerow(headerRow)

   for fn in fns:
      #basename = os.path.basename(fn)
      basename, ext = os.path.splitext(os.path.basename(fn))
      if basename.find('Neg') > 0:
         mult = -1.0
      else:
         mult = 1.0
      moransI = mult * float(basename[6:])
      
      asciiSize = getFileSize(fn)
      
      tiffSize = getFileSize(getFilename(DATA_DIR, TIFF, basename, 'tif'))
      
      # Binary
      binFn = getFilename(DATA_DIR, BIN, basename, 'bin')
      binMtx = convertAscToMatrix(fn, threshold=2)
      binMtx.writeFile(binFn)
      binSize = getFileSize(binFn)
      
      # Num presences & variance
      numPresences, variance = getNumPresencesAndVariance(binMtx.data)
      
      # RLE
      rleFn = getFilename(DATA_DIR, RLE, basename, 'rle')
      rleMtx = NormalRLECompressedMatrix()
      rleMtx.compress(binMtx)
      rleMtx.writeFile(rleFn)
      rleSize = getFileSize(rleFn)
      
      # Hilbert
      hilbFn = getFilename(DATA_DIR, HILBERT, basename, 'hilb')
      hilbMtx = HilbertRLECompressedMatrix()
      hilbMtx.compress(binMtx)
      hilbMtx.writeFile(hilbFn)
      hilbertSize = getFileSize(hilbFn)
   
      # QTree
      qTreeFn = getFilename(DATA_DIR, QTREE, basename, 'qtree')
      qMtx = QuadtreeCompressedMatrix()
      qMtx.compress(binMtx)
      qMtx.writeFile(qTreeFn)
      qTreeSize = getFileSize(qTreeFn)
      
      # STree
      sTreeFn = getFilename(DATA_DIR, STREE, basename, 'stree')
      sMtx = STreeCompressedMatrix()
      sMtx.compress(binMtx)
      sMtx.writeFile(sTreeFn)
      sTreeSize = getFileSize(sTreeFn)
   
      # Combo 6
      c6Fn = getFilename(DATA_DIR, COMBO6, basename, 'combo6')
      c6Mtx = QuadtreeHilbertCompressedMatrix(rleThreshold=6)
      c6Mtx.compress(binMtx)
      c6Mtx.writeFile(c6Fn)
      combo6Size = getFileSize(c6Fn)
      
      # Combo 8
      c8Fn = getFilename(DATA_DIR, COMBO8, basename, 'combo8')
      c8Mtx = QuadtreeHilbertCompressedMatrix(rleThreshold=8)
      c8Mtx.compress(binMtx)
      c8Mtx.writeFile(c8Fn)
      combo8Size = getFileSize(c8Fn)
   
      # Combo 10
      c10Fn = getFilename(DATA_DIR, COMBO10, basename, 'combo10')
      c10Mtx = QuadtreeHilbertCompressedMatrix(rleThreshold=10)
      c10Mtx.compress(binMtx)
      c10Mtx.writeFile(c10Fn)
      combo10Size = getFileSize(c10Fn)

      row = [ moransI, numPresences, variance, asciiSize, tiffSize, binSize, 
             rleSize, hilbertSize, qTreeSize, sTreeSize, combo6Size, 
             combo8Size, combo10Size]
      writer.writerow(row)
