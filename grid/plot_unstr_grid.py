import numpy as np  
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr

def read_gmsh(filename):
    
    f = open(filename,'r')
    lines = f.read().splitlines()

    lines.pop(0)
    lines.pop(0)
    lines.pop(0)
    lines.pop(0)
    nn = int(lines.pop(0))

    # Initialize node coordinate and depth arrays
    xy = np.zeros((nn,2),dtype=np.double)
    depth = np.zeros((nn,),dtype=np.double)
  
    # Read coordinates and depths
    for i in range(nn):
      line = lines.pop(0).split()
      j = int(line[0])-1
      xy[j,0] = float(line[1])
      xy[j,1] = float(line[2])
      depth[j] = float(line[3])

    lines.pop(0)
    lines.pop(0)
    ne = int(lines.pop(0))

    ect = np.zeros((ne,3),dtype=np.int32)

    for i in range(ne): 

      line = lines.pop(0).split()
      j = int(line[0])-1
      ect[j,0] = int(line[6])-1
      ect[j,1] = int(line[7])-1
      ect[j,2] = int(line[8])-1

    return xy,ect

def plot_field(xy,ect,ds):

  fig = plt.figure(figsize=[18.0,9.0])
  #ax = fig.add_subplot(1,1,1, projection=ccrs.PlateCarree())
  ax = fig.add_subplot(1,1,1, projection=ccrs.NorthPolarStereo())
  ax.set_extent([-180.0,180.0,85.0,90.0], crs=ccrs.PlateCarree())


  cf = ax.triplot(xy[:,0],xy[:,1],ect,  transform=ccrs.Geodetic())
  ax.scatter(ds['grid_center_lon'].values,ds['grid_center_lat'].values,marker='.',color='grey', transform=ccrs.PlateCarree())
  ax.scatter(ds['grid_corner_lon'].values,ds['grid_corner_lat'].values,marker='.',color='k', transform=ccrs.PlateCarree())
  ax.scatter(0.0,90.0,marker='.',color='r',transform=ccrs.PlateCarree())
  ax.add_feature(cfeature.LAND, alpha=0.5, zorder=100)
  ax.add_feature(cfeature.LAKES, alpha=0.5, zorder=101)
  ax.add_feature(cfeature.COASTLINE, zorder=101)
  plt.savefig('test.png')
  plt.close()



if __name__ == "__main__":

  xy,ect = read_gmsh('mesh.msh')
  ds = xr.open_dataset('scrip.nc')
  print(ds['grid_center_lon'].values)
  print(ds['grid_center_lat'].values)

  ds['grid_center_lon'].values = np.where(ds['grid_center_lon'].values > 180.0, ds['grid_center_lon'].values-360, ds['grid_center_lon'].values)
  print(ds['grid_center_lon'].values)
  print(ds['grid_center_lat'].values)

  plot_field(xy,ect,ds)
