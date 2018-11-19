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

##################################################################################################
##################################################################################################

def interpolate_data_to_grid(grid_file,data_file,var):

  # Open files  
  data_nc = netCDF4.Dataset(data_file,'r')
  grid_nc = netCDF4.Dataset(grid_file,'r')

  # Get grid from data file
  lon_data = data_nc.variables['lon'][:]
  lat_data = data_nc.variables['lat'][:]
  time = data_nc.variables['time'][:]
  nsnaps = time.size
  nlon = lon_data.size
  nlat = lat_data.size
  data = np.zeros((nsnaps,nlat,nlon))

  # Get grid from grid file
  lon_grid = grid_nc.variables['lon'][:]
  lat_grid = grid_nc.variables['lat'][:]
  Lon,Lat = np.meshgrid(lon_grid,lat_grid)
  grid_points = np.column_stack((Lon.ravel(),Lat.ravel()))
  nlon = lon_grid.size
  nlat = lat_grid.size
  interp_data = np.zeros((nsnaps,nlat,nlon))
  print interp_data.shape

  for i,t in enumerate(time):
    print i

    data[i,:,:] = data_nc.variables[var][i,:,:]

    interpolator = interpolate.RegularGridInterpolator((lon_data,lat_data),data[i,:,:].T)
    interp_data[i,:,:] = np.reshape(interpolator(grid_points),(nlat,nlon))

  return lon_grid,lat_grid,interp_data,lon_data,lat_data,data

##################################################################################################
##################################################################################################

def plot_interp_data(lon_data,lat_data,data,lon_grid,lat_grid,interp_data,var_label):

  fig = plt.figure()
  ax0 = fig.add_subplot(2,1,1)
  levels = np.linspace(np.amin(data),np.amax(data),10)
  cf = ax0.contourf(lon_data,lat_data,data,levels=levels)
  ax0.set_title('data')
  cbar = fig.colorbar(cf,ax=ax0)
  cbar.set_label(var_label)
  ax1 = fig.add_subplot(2,1,2)
  cf = ax1.contourf(lon_grid,lat_grid,interp_data,levels=levels)
  ax1.set_title('interpolated data')
  cbar = fig.colorbar(cf,ax=ax1)
  cbar.set_label(var_label)
  fig.tight_layout()
  fig.savefig('vel_'+str(i).zfill(3)+'.png',box_inches='tight')
  plt.close()

##################################################################################################
##################################################################################################

if __name__ == '__main__':

  grid_file = '/users/sbrus/scratch3/ACME/input_data/ocn/iaf/ncep.u_10.T62.1948.nc'
  data_file = '/users/sbrus/scratch4/WW3_testing/wind_data/june/wnd10mx0.5.gdas.200506.ww3.nc'

  lon_grid,lat_grid,u_interp,lon_data,lat_data,u_data = interpolate_data_to_grid(grid_file,data_file,'U_GRD_L103')
  lon_grid,lat_grid,v_interp,lon_data,lat_data,v_data = interpolate_data_to_grid(grid_file,data_file,'V_GRD_L103')
 
 
  for i in range(u_data.shape[0]):
    print i 

    data = np.sqrt(np.square(u_data[i,:,:]) + np.square(v_data[i,:,:]))
    interp_data = np.sqrt(np.square(u_interp[i,:,:]) + np.square(v_interp[i,:,:]))
    
    plot_interp_data(lon_data,lat_data,data,lon_grid,lat_grid,interp_data,'velocity magnitude')
