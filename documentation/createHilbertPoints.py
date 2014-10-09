from osgeo import ogr

myWkt = [
'POINT(5 5)',
'POINT(5 15)',
'POINT(15 15)',
'POINT(15 5)',
'POINT(25 5)',
'POINT(35 5)',
'POINT(35 15)',
'POINT(25 15)',
'POINT(25 25)',
'POINT(35 25)',
'POINT(35 35)',
'POINT(25 35)',
'POINT(15 35)',
'POINT(15 25)',
'POINT(5 25)',
'POINT(5 35)',
'POINT(5 45)',
'POINT(15 45)',
'POINT(15 55)',
'POINT(5 55)',
'POINT(5 65)',
'POINT(5 75)',
'POINT(15 75)',
'POINT(15 65)',
'POINT(25 65)',
'POINT(25 75)',
'POINT(35 75)',
'POINT(35 65)',
'POINT(35 55)',
'POINT(25 55)',
'POINT(25 45)',
'POINT(35 45)',
'POINT(45 45)',
'POINT(55 45)',
'POINT(55 55)',
'POINT(45 55)',
'POINT(45 65)',
'POINT(45 75)',
'POINT(55 75)',
'POINT(55 65)',
'POINT(65 65)',
'POINT(65 75)',
'POINT(75 75)',
'POINT(75 65)',
'POINT(75 55)',
'POINT(65 55)',
'POINT(65 45)',
'POINT(75 45)',
'POINT(75 35)',
'POINT(75 25)',
'POINT(65 25)',
'POINT(65 35)',
'POINT(55 35)',
'POINT(45 35)',
'POINT(45 25)',
'POINT(55 25)',
'POINT(55 15)',
'POINT(45 15)',
'POINT(45 5)',
'POINT(55 5)',
'POINT(65 5)',
'POINT(65 15)',
'POINT(75 15)',
'POINT(75 5)',
]


if __name__ == "__main__":
   driverName = "ESRI Shapefile"
   fn = "/home/cjgrady/Desktop/hilbertPoints.shp"

   driver = ogr.GetDriverByName(driverName)
   ds = driver.CreateDataSource(fn)
   layer = ds.CreateLayer('hilbertPoints', geom_type=ogr.wkbPoint)
   layer.CreateField(ogr.FieldDefn('index', ogr.OFTInteger))
   layerDfn = layer.GetLayerDefn()

   i = 1
   for txt in myWkt:
      point = ogr.CreateGeometryFromWkt(txt)

      feature = ogr.Feature(layerDfn)
      feature.SetGeometry(point)
      feature.SetFID(i)
      feature.SetField('index', i)
      i = i+1
      layer.CreateFeature(feature)

   ds.Destroy()
