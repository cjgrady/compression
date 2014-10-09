from osgeo import ogr

data = [
        [1, 1, 1, 1, 1, 1, 1, 1],
        [2, 2, 1, 1, 2, 2, 1, 1],
        [3, 2, 1, 1, 2, 1, 1, 2],
        [3, 2, 2, 2, 2, 1, 1, 4],
        [3, 3, 3, 1, 1, 1, 2, 4],
        [3, 3, 3, 1, 1, 2, 2, 4],
        [3, 3, 3, 4, 4, 4, 2, 4],
        [3, 3, 3, 3, 4, 4, 4, 4]
       ]

if __name__ == "__main__":
   driverName = "ESRI Shapefile"
   fn = "/home/cjgrady/Desktop/hilbertCells.shp"

   driver = ogr.GetDriverByName(driverName)
   ds = driver.CreateDataSource(fn)
   layer = ds.CreateLayer('hilbertCells', geom_type=ogr.wkbPolygon)
   layer.CreateField(ogr.FieldDefn('category', ogr.OFTInteger))
   layerDfn = layer.GetLayerDefn()

   i = 1
   cellSize = 10
   for y in xrange(len(data)):
      for x in xrange(len(data[y])):
         x0 = x*cellSize
         x1 = (x+1) * cellSize
         y0 = (7-y) * cellSize
         y1 = (8-y) * cellSize
         wkt = "POLYGON(({x0} {y0}, {x1} {y0}, {x1} {y1}, {x0} {y1}, {x0} {y0}))".format(x0=x0, x1=x1, y0=y0, y1=y1)
         print wkt
         
         poly = ogr.CreateGeometryFromWkt(wkt)
         print " created polygon geometry"

         feature = ogr.Feature(layerDfn)
         print "feature"
         feature.SetGeometry(poly)
         print "set geometry"
         feature.SetFID(i)
         print "set fid"
         feature.SetField('category', data[y][x])
         print "set field category"
         i = i+1
         layer.CreateFeature(feature)
         print "created feature"

   ds.Destroy()
