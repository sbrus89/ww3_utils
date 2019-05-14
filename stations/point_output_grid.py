import numpy as np

ddeg = 3.0

lon_min = -180.0
lon_max =  180.0
lon = np.arange(lon_min,lon_max,ddeg)
print lon
print lon.size

lat_min = -80.0
lat_max =  80.0
lat = np.arange(lat_min,lat_max,ddeg)
print lat
print lat.size

lon_grid,lat_grid = np.meshgrid(lon,lat)

lon_pts = np.ravel(lon_grid)
lat_pts = np.ravel(lat_grid)

f = open('./stations.txt','a')
for i in range(lon_pts.size):
  lon = lon_pts[i]
  lat = lat_pts[i]

  name = str(lon)+','+str(lat)

  print name
  f.write(str(lon)+" "+str(lat)+" '"+name+"'\n")

f.close()  

