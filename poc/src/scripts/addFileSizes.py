"""
@summary: Adds file sizes to csv file
"""
import csv
import os

inCsvFn = '/data/geo716/final project/data/projections2.csv'
outCsvFn = '/data/geo716/final project/data/projections.csv'

binDir = '/data/geo716/final project/data/bins/'
rleDir = '/data/geo716/final project/data/rles/'
hilbertDir = '/data/geo716/final project/data/hilberts/'


with open(inCsvFn) as csvIn:
   with open(outCsvFn, 'w') as csvOut:
      reader = csv.reader(csvIn)
      writer = csv.writer(csvOut)
      
      headers = reader.next()
      headers.append('Binary size')
      headers.append('RLE size')
      headers.append('Hilbert size')

      all = []
      all.append(headers)
      
      i = 0
      for row in reader:
         binFn = os.path.join(binDir, '%s.bin' % i)
         rleFn = os.path.join(rleDir, '%s.rle' % i)
         hilbertFn = os.path.join(hilbertDir, '%s.hilb' % i)
         
         row.append(os.path.getsize(binFn))
         row.append(os.path.getsize(rleFn))
         row.append(os.path.getsize(hilbertFn))

         all.append(row)         
         i += 1
      writer.writerows(all)
      
