from netCDF4 import Dataset
import numpy as np
import shapely
from scipy import interpolate
from geometric_features import read_feature_collection
import calc_distance

def create_initial_points(meshfile,lon,lat,hfunction,outfile):

  # Open MPAS mesh and get cell variables
  nc_file = Dataset(meshfile,'r')
  lonCell = nc_file.variables['lonCell'][:]
  latCell = nc_file.variables['latCell'][:]
  bottomDepth = nc_file.variables['bottomDepth'][:]

  # Transform 0,360 range to -180,180 
  idx, = np.where(lonCell > np.pi) 
  lonCell[idx] = lonCell[idx] - 2.0*np.pi

  # Interpolate hfunction onto mesh cell centers
  hfun = interpolate.RegularGridInterpolator((np.radians(lon),np.radians(lat)),hfunction.T)
  mesh_pts = np.vstack((lonCell,latCell)).T
  hfun_interp = hfun(mesh_pts)

  # Find cells in refined region of waves mesh
  max_res = np.amax(hfunction)
  idx, = np.where(hfun_interp < 0.5*max_res )

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
