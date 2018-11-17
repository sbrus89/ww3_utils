import netCDF4 
import matplotlib.pyplot as plt
import numpy as np
import glob
import pprint
import datetime
import os
import yaml
import subprocess
from mpl_toolkits.basemap import Basemap
from scipy import interpolate
plt.switch_backend('agg')


def interpolate_data_to_grid(grid_file,data_file,var):
  
  data_nc = netCDF4.Dataset(data_file,'r')
  grid_nc = netCDF4.Dataset(grid_file,'r')

  lon_data = data_nc.variables['lon'][:]
  lat_data = data_nc.variables['lat'][:]
  time = data_nc.variables['time'][:]
  nsnaps = time.size

  lon_grid = grid_nc.variables['lon'][:]
  lat_grid = grid_nc.variables['lat'][:]
  Lon,Lat = np.meshgrid(lon_grid,lat_grid)
  grid_points = np.column_stack((Lon.ravel(),Lat.ravel()))
  nlon = lon_grid.size
  nlat = lat_grid.size

  interp_data = np.zeros((nsnaps,nlon,nlat))
  print interp_data.shape

  for i,t in enumerate(time):
    print i
    data = data_nc.variables[var][i,:,:]

    interpolator = interpolate.RegularGridInterpolator((lon_data,lat_data),data.T)
    interp_data[i,:,:] = np.reshape(interpolator(grid_points),(nlon,nlat))

  return interp_data

if __name__ == '__main__':

  grid_file = '/users/sbrus/scratch3/ACME/input_data/ocn/iaf/ncep.u_10.T62.1948.nc'
  data_file = '/users/sbrus/scratch4/WW3_testing/wind_data/june/wnd10mx0.5.gdas.200506.ww3.nc'

  u_interp = interpolate_data_to_grid(grid_file,data_file,'U_GRD_L103')
  v_interp = interpolate_data_to_grid(grid_file,data_file,'V_GRD_L103')
 

