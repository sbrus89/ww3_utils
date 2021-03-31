from netCDF4 import Dataset
import numpy as np
from scipy import interpolate
import calc_distance

km = 1000.0

def create_initial_points(meshfile,lon,lat,hfunction,sphere_radius,outfile):

  # Open MPAS mesh and get cell variables
  nc_file = Dataset(meshfile,'r')
  lonCell = nc_file.variables['lonCell'][:]
  latCell = nc_file.variables['latCell'][:]

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

  # Find boundary cells
  nEdgesOnCell = nc_file.variables['nEdgesOnCell'][:]
  cellsOnCell = nc_file.variables['cellsOnCell'][:]
  nCellsOnCell = np.count_nonzero(cellsOnCell,axis=1)
  is_boundary_cell = np.equal(nCellsOnCell,nEdgesOnCell)
  idx_bnd, = np.where(is_boundary_cell == False)

  # Force inclusion of all boundary cells
  idx = np.union1d(idx,idx_bnd)

  # Get coordinates of points
  lon = lonCell[idx]
  lat = latCell[idx]
  npt = idx.size

  # Change to Cartesian coordinates
  x,y,z = calc_distance.lonlat2xyz(lon,lat,sphere_radius)

  # Get coordinates and ID into structured array (for use with np.savetxt)
  pt_list = []
  for i in range(npt):
    pt_list.append((x[i],y[i],z[i],-1))  # ID of -1 specifies that node is fixed
  pt_type = np.dtype({'names':['x','y','z','id'],'formats':[np.float64, np.float64, np.float64, np.int32]})
  pts = np.array(pt_list,dtype=pt_type)

  # Write initial conditions file
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
