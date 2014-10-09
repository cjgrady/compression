"""
@summary: Get the number of presences in a layer
"""
import csv
import os
from matrix.matrix import Matrix

inDir = '/data/geo716/final project/data/bins'
inCsvFn = '/data/geo716/final project/data/projectionsNew.csv'
outCsvFn = '/data/geo716/final project/data/projections2.csv'

with open(inCsvFn) as csvIn:
   with open(outCsvFn, 'w') as csvOut:
      reader = csv.reader(csvIn)
      writer = csv.writer(csvOut)
      
      headers = reader.next()
      headers.append('Number of presence')

      all = []
      all.append(headers)
      
      i = 0
      for row in reader:
         fn = os.path.join(inDir, '%s.bin' % i)
         mtx = Matrix()
         mtx.readFile(fn)
         row.append(sum([sum(r) for r in mtx.data]))
         all.append(row)         
         i += 1
      writer.writerows(all)
      
