from netCDF4 import Dataset
import numpy as np
import shapely
from geometric_features import read_feature_collection
import calc_distance

def create_initial_points(mesh,outfile):

  # Open MPAS mesh and get cell variables
  nc_file = Dataset(mesh,'r')
  lonCell = nc_file.variables['lonCell'][:]
  latCell = nc_file.variables['latCell'][:]
  bottomDepth = nc_file.variables['bottomDepth'][:]

  # Transform 0,360 range to -180,180 
  idx, = np.where(lonCell > np.pi) 
  lonCell[idx] = lonCell[idx] - 2.0*np.pi

  ## Find points inside geojson regions
  #fileName = 'map'
  #in_points = []
  #fc = read_feature_collection('{}.geojson'.format(fileName))
  #for feature in fc.features:
  #  shape = shapely.geometry.Polygon(feature['geometry']['coordinates'][0]) 
  #  for i in range(lonCell.size):
  #    point = shapely.geometry.Point(np.degrees(lonCell[i]),np.degrees(latCell[i]))
  #    test = shape.contains(point)
  #    if (test == True) and (i not in in_points):
  #      print(i)
  #      in_points.append(i)
  #in_points = np.array(in_points)

  shpfiles = [ 
               "/home/sbrus/run/WW3_unstructured/OceanMesh2D/utilities/GSHHS/c/GSHHS_c_L1.shp",
               "/home/sbrus/run/WW3_unstructured/OceanMesh2D/utilities/GSHHS/c/GSHHS_c_L6.shp"
           ]
  D = calc_distance.distance_to_shapefile_points(shpfiles,lonCell,latCell)

  ## Find points shallower then depth threshold
  #above_depth, = np.where(bottomDepth < 4000.0)
  ## Find points that are inside regions and shallower than threshold
  #idx = np.intersect1d(above_depth,in_points)

  idx, = np.where((bottomDepth < 4000.0) & (D < 1000.0))
  #idx, = np.where(D < 1000.0)

  # Get coordinates of points
  lon = lonCell[idx]
  lat = latCell[idx]
  npt = idx.size

  # Change to Cartesian coordinates
  x,y,z = calc_distance.lonlat2xyz(lon,lat)

  # Specify that initial points are fixed
  ID = -1*np.ones(x.shape)

  pt_list = []
  for i in range(npt):
    pt_list.append((x[i]/1000.0,y[i]/1000.0,z[i]/1000.0,-1))

  pt_type = np.dtype({'names':['x','y','z','id'],'formats':[np.float64, np.float64, np.float64, np.int32]})
  pts = np.array(pt_list,dtype=pt_type)

  f = open(outfile,'w')
  f.write('# Initial coordinates \n')
  f.write('MSHID=3;EUCLIDEAN-MESH\n')
  f.write('NDIMS=3\n')
  f.write('POINT='+str(npt)+'\n')
  np.savetxt(f,pts,fmt='%.12e;%.12e;%.12e;%2i')
  f.close()

  return 


if __name__ == '__main__':
  create_initial_points('ocean.WC14to60E2r3.200714_scaled.nc','init.msh')
