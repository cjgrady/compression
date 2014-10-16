"""
@summary: This script will read snow cover data and compress it
@author: CJ Grady
@version: 1.0
@status: alpha

@license: gpl2
@copyright: Copyright (C) 2014, University of Kansas Center for Research

          Lifemapper Project, lifemapper [at] ku [dot] edu, 
          Biodiversity Institute,
          1345 Jayhawk Boulevard, Lawrence, Kansas, 66045, USA
   
          This program is free software; you can redistribute it and/or modify 
          it under the terms of the GNU General Public License as published by 
          the Free Software Foundation; either version 2 of the License, or (at 
          your option) any later version.
  
          This program is distributed in the hope that it will be useful, but 
          WITHOUT ANY WARRANTY; without even the implied warranty of 
          MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
          General Public License for more details.
  
          You should have received a copy of the GNU General Public License 
          along with this program; if not, write to the Free Software 
          Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
          02110-1301, USA.
"""
import os
from matrix.matrix import Grid
from methods.combo.quadTreeCombo import QuadtreeComboCompressedGrid

from methods.rle.hilbert import HilbertRLECompressedGrid
from methods.rle.normal import NormalRLECompressedGrid

INPUT_FN = '../../data/snow2.txt'
HILB_OUTPUT_FN = '../../data/snow2-2.hilb'
NORM_OUTPUT_FN = '../../data/snow2-2.rle'
OUT_DIR = '../../data/'

# .............................................................................
if __name__ == "__main__":
   
   if os.path.exists(INPUT_FN):
      # Read in snow data
      data = []
      with open(INPUT_FN) as inF:
         for line in inF:
            data.append([int(i) for i in line.split(',')])

      # Convert to a grid
      grid = Grid(data)
      
      # Feed into compression
      cmp1 = NormalRLECompressedGrid(grid=grid)
      cmp2 = HilbertRLECompressedGrid(grid=grid)
      cmp3 = QuadtreeComboCompressedGrid(NormalRLECompressedGrid, orderThreshold=6, grid=grid)
      cmp4 = QuadtreeComboCompressedGrid(HilbertRLECompressedGrid, orderThreshold=6, grid=grid)
      cmp5 = QuadtreeComboCompressedGrid(NormalRLECompressedGrid, orderThreshold=8, grid=grid)
      cmp6 = QuadtreeComboCompressedGrid(HilbertRLECompressedGrid, orderThreshold=8, grid=grid)
      cmp7 = QuadtreeComboCompressedGrid(NormalRLECompressedGrid, orderThreshold=10, grid=grid)
      cmp8 = QuadtreeComboCompressedGrid(HilbertRLECompressedGrid, orderThreshold=10, grid=grid)
      cmp9 = QuadtreeComboCompressedGrid(NormalRLECompressedGrid, orderThreshold=2, grid=grid)
      cmp10 = QuadtreeComboCompressedGrid(HilbertRLECompressedGrid, orderThreshold=2, grid=grid)
      cmp11 = QuadtreeComboCompressedGrid(NormalRLECompressedGrid, orderThreshold=4, grid=grid)
      cmp12 = QuadtreeComboCompressedGrid(HilbertRLECompressedGrid, orderThreshold=4, grid=grid)
      cmp13 = QuadtreeComboCompressedGrid(NormalRLECompressedGrid, orderThreshold=12, grid=grid)
      cmp14 = QuadtreeComboCompressedGrid(HilbertRLECompressedGrid, orderThreshold=12, grid=grid)
      cmp15 = QuadtreeComboCompressedGrid(NormalRLECompressedGrid, orderThreshold=14, grid=grid)
      cmp16 = QuadtreeComboCompressedGrid(HilbertRLECompressedGrid, orderThreshold=14, grid=grid)
      
      # Write out compressed data      
      cmp1.write(NORM_OUTPUT_FN)
      cmp2.write(HILB_OUTPUT_FN)
      
      # Combo tests
      cmp3.write(os.path.join(OUT_DIR, 'snow2-combo-normal-6.bin'))
      cmp4.write(os.path.join(OUT_DIR, 'snow2-combo-hilb-6.bin'))
      cmp5.write(os.path.join(OUT_DIR, 'snow2-combo-normal-8.bin'))
      cmp6.write(os.path.join(OUT_DIR, 'snow2-combo-hilb-8.bin'))
      cmp7.write(os.path.join(OUT_DIR, 'snow2-combo-normal-10.bin'))
      cmp8.write(os.path.join(OUT_DIR, 'snow2-combo-hilb-10.bin'))
      cmp9.write(os.path.join(OUT_DIR, 'snow2-combo-normal-2.bin'))
      cmp10.write(os.path.join(OUT_DIR, 'snow2-combo-hilb-2.bin'))
      cmp11.write(os.path.join(OUT_DIR, 'snow2-combo-normal-4.bin'))
      cmp12.write(os.path.join(OUT_DIR, 'snow2-combo-hilb-4.bin'))
      cmp13.write(os.path.join(OUT_DIR, 'snow2-combo-normal-12.bin'))
      cmp14.write(os.path.join(OUT_DIR, 'snow2-combo-hilb-12.bin'))
      cmp15.write(os.path.join(OUT_DIR, 'snow2-combo-normal-14.bin'))
      cmp16.write(os.path.join(OUT_DIR, 'snow2-combo-hilb-14.bin'))
      
   else:
      print "File %s does not exist" % INPUT_FN
