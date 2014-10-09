"""
@summary: Creates a CSV file with all of the measurements and stats collected
@author: CJ Grady
"""
import csv
import os

from matrix.matrix import Matrix
from tools.moransI import moransI

BASE_DIR = "/data/geo716/final project/data"
BASE_CSV = "/data/geo716/final project/data/base.csv"
OUT_CSV = "/data/geo716/final project/data/out.csv"

# .............................................................................
def calculateMoransI(mtx):
   c = [
        [0, 1, 0],
        [1, 0, 1],
        [0, 1, 0]
       ]
   return moransI(mtx, c)

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
   
# .............................................................................
def getFileName(directory, name, extension):
   return os.path.join(BASE_DIR, directory, '%s.%s' % (name, extension))

# .............................................................................
def getFileSize(fn):
   return os.path.getsize(fn)

# .............................................................................
if __name__ == "__main__":
   # Read base CSV file
   values = []
   with open(BASE_CSV) as csvInFile:
      reader = csv.reader(csvInFile)
      for row in reader:
         values.append(row)

   with open(OUT_CSV, 'w') as csvOutFile:
      writer = csv.writer(csvOutFile)
      headerRow = [
                   'Id', 'Projection Id', 'Occurrence Set Id', 'Scenario Id',
                   'Algorithm Code', 'Number of Presences', 'Morans I',
                   'Variance', 'ASCII Size', 'Tiff Size', 'Binary Size',
                   'RLE Size', 'Hilbert Size', 'Quadtree Size', 'S-Tree Size',
                   'Combo 6 Size', 'Combo 8 Size', 'Combo 10 Size'
                  ]
      writer.writerow(headerRow)
      for idStr, prjId, occId, scnId, algoCode, tiffSize in values[1:]:
         print idStr
         id = int(idStr)
         mtx = Matrix()
         mtx.readFile(getFileName('bins', id, 'bin'))
         numPres, variance = getNumPresencesAndVariance(mtx.data)
         mI = calculateMoransI(mtx.data)
         asciiSize = getFileSize(getFileName('ascs', id, 'asc'))
         tiffSize = tiffSize
         binSize = getFileSize(getFileName('bins', id, 'bin'))
         rleSize = getFileSize(getFileName('rles', id, 'rle'))
         hilbertSize = getFileSize(getFileName('hilberts', id, 'hilb'))
         qTreeSize = getFileSize(getFileName('qtrees', id, 'qtree'))
         sTreeSize = getFileSize(getFileName('strees', id, 'stree'))
         combo6Size = getFileSize(getFileName('combos6', id, 'combo6'))
         combo8Size = getFileSize(getFileName('combos8', id, 'combo8'))
         combo10Size = getFileSize(getFileName('combos10', id, 'combo10'))
         
         row = [id, prjId, occId, scnId, algoCode, numPres, mI, variance,
                asciiSize, tiffSize, binSize, rleSize, hilbertSize,
                qTreeSize, sTreeSize, combo6Size, combo8Size, combo10Size]
         writer.writerow(row)
