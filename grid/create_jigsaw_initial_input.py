from netCDF4 import Dataset
import numpy as np

def create_initial_points(mesh,outfile):

  nc_file = Dataset(mesh,'r')
  lonCell = nc_file.variables['lonCell'][:]
  latCell = nc_file.variables['latCell'][:]
  bottomDepth = nc_file.variables['bottomDepth'][:]

  idx, = np.where(bottomDepth < 4000.0)
 
  pt_list = [] 
  for i in idx:
    lon = lonCell[i]
    lat = latCell[i]

    if lon > np.pi:
      lon = lon - 2.0*np.pi

    R = 6371.0
    x = R*np.cos(lat)*np.cos(lon) 
    y = R*np.cos(lat)*np.sin(lon)
    z = R*np.sin(lat)  

    pt_list.append((x,y,z,-1))
    #pt_list.append((lon,lat,-1))


  pt_type = np.dtype({'names':['x','y','z','id'],'formats':[np.float64, np.float64, np.float64, np.int32]})
  #pt_type = np.dtype({'names':['lon','lat','id'],'formats':[np.float64, np.float64, np.int32]})
  pts = np.array(pt_list,dtype=pt_type)

  f = open(outfile,'w')
  f.write('# Initial coordinates \n')
  f.write('MSHID=3;EUCLIDEAN-MESH\n')
  #f.write('MSHID=3;ELLIPSOID-MESH\n')
  #f.write('RADII=6371;6371;6371\n')
  f.write('NDIMS=3\n')
  #f.write('NDIMS=2\n')
  npt = len(pt_list)
  f.write('POINT='+str(npt)+'\n')
  np.savetxt(f,pts,fmt='%.12e;%.12e;%.12e;%2i')
  #np.savetxt(f,pts,fmt='%.12e;%.12e;%2i')
  f.close()

  return 


if __name__ == '__main__':
  create_initial_points('ocean.WC14to60E2r3.200714_scaled.nc','init.msh')
