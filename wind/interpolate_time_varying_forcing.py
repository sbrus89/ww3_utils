import netCDF4 
import matplotlib.pyplot as plt
import numpy as np
import glob
import pprint
import datetime
import os
import yaml
import subprocess
import argparse
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
  lon_data = np.append(lon_data,360.0)
  lat_data = np.flipud(data_nc.variables['lat'][:])
  time = data_nc.variables['time'][:]
  nsnaps = time.size
  nlon = lon_data.size
  nlat = lat_data.size
  data = np.zeros((nsnaps,nlat,nlon))
  print data.shape

  # Get grid from grid file
  lon_grid = grid_nc.variables['lonCell'][:]*180.0/np.pi
  lat_grid = grid_nc.variables['latCell'][:]*180.0/np.pi
  grid_points = np.column_stack((lon_grid,lat_grid))
  ncells = lon_grid.size
  interp_data = np.zeros((nsnaps,ncells))
  print interp_data.shape
  print np.amin(lon_grid),np.amax(lon_grid)
  print np.amin(lat_grid),np.amax(lat_grid)

  # Interpolate timesnaps
  for i,t in enumerate(time):
    print i

    # Get data to interpolate
    data[i,:,0:-1] = np.flipud(data_nc.variables[var][i,:,:])
    data[i,:,-1] = data[i,:,0]

    # Interpolate data onto new grid
    interpolator = interpolate.RegularGridInterpolator((lon_data,lat_data),data[i,:,:].T,bounds_error=False,fill_value=0.0)
    interp_data[i,:] = interpolator(grid_points)

  return lon_grid,lat_grid,interp_data,lon_data,lat_data,data

##################################################################################################
##################################################################################################

def write_to_file(filename,data,var):

  data_nc = netCDF4.Dataset(filename,'a')

  # Find dimesions
  ncells = data.shape[1]
  nsnaps = data.shape[0]

  # Declare dimensions
  data_nc.createDimension('nCells',ncells)
  data_nc.createDimension('StrLen',64)
  data_nc.createDimension('Time',None)

  # Declear variables
  data_var = data_nc.createVariable(var,np.float64,('Time','nCells'))
     #time = data_nc.createVariable('xtime',np.float64,('Time','StrLen'))

  # Set variables
  data_var[:,:] = data
  data_nc.close()

##################################################################################################
##################################################################################################

def plot_interp_data(lon_data,lat_data,data,lon_grid,lat_grid,interp_data,var_label):
  print i

  levels = np.linspace(np.amin(data),np.amax(data),10)

  fig = plt.figure()
  ax0 = fig.add_subplot(2,1,1)
  cf = ax0.contourf(lon_data,lat_data,data,levels=levels)
  ax0.set_title('data')
  cbar = fig.colorbar(cf,ax=ax0)
  cbar.set_label(var_label)

  ax1 = fig.add_subplot(2,1,2)
  cf = ax1.tricontourf(lon_grid,lat_grid,interp_data,levels=levels)
  ax1.set_title('interpolated data')
  cbar = fig.colorbar(cf,ax=ax1)
  cbar.set_label(var_label)

  fig.tight_layout()
  fig.savefig('vel_'+str(i).zfill(3)+'.png',box_inches='tight')
  plt.close()

##################################################################################################
##################################################################################################

if __name__ == '__main__':

  parser = argparse.ArgumentParser()
  parser.add_argument('--plot',action='store_true')
  args = parser.parse_args()
  
  # Files to interpolate to/from
  grid_file = '/users/sbrus/scratch4/MPAS-O_testing/ocean/global_ocean/USDEQU120cr10rr1/init/culled_mesh/culled_mesh.nc'
  data_file = '/users/sbrus/scratch4/MPAS-O_testing/time_varying_forcing/wind_data/wnd10m.cdas1.201210.grb2.nc'

  # Interpolation of u and v velocities
  lon_grid,lat_grid,u_interp,lon_data,lat_data,u_data = interpolate_data_to_grid(grid_file,data_file,'U_GRD_L103')
  lon_grid,lat_grid,v_interp,lon_data,lat_data,v_data = interpolate_data_to_grid(grid_file,data_file,'V_GRD_L103')
 
  # Calculate and plot velocity magnitude
  if args.plot:
    for i in range(u_data.shape[0]):
      print i 
  
      data = np.sqrt(np.square(u_data[i,:,:]) + np.square(v_data[i,:,:]))
      interp_data = np.sqrt(np.square(u_interp[i,:]) + np.square(v_interp[i,:]))
      
      plot_interp_data(lon_data,lat_data,data,lon_grid,lat_grid,interp_data,'velocity magnitude')

  write_to_file('test.nc',u_interp,'u_10')
  write_to_file('test.nc',v_interp,'v_10')
