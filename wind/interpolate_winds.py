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

  # Interpolate timesnaps
  for i,t in enumerate(time):
    print i

    # Get data to interpolate
    data[i,:,:] = data_nc.variables[var][i,:,:]

    # Interpolate data onto new grid
    interpolator = interpolate.RegularGridInterpolator((lon_data,lat_data),data[i,:,:].T)
    interp_data[i,:,:] = np.reshape(interpolator(grid_points),(nlat,nlon))

  return lon_grid,lat_grid,interp_data,lon_data,lat_data,data

##################################################################################################
##################################################################################################

def write_to_file(filename,data,var):

  data_nc = netCDF4.Dataset(filename,'w')

  # Find dimesions
  nlat = data.shape[1]
  nlon = data.shape[2]
  nsnaps = data.shape[0]

  # Declare dimensions
  data_nc.createDimension('lat',nlat)
  data_nc.createDimension('lon',nlon)
  data_nc.createDimension('time',None)

  # Declear variables
  data_var = data_nc.createVariable(var,np.float64,('time','lat','lon'))
  lat = data_nc.createVariable('lat',np.float64,('lat',))
  lon = data_nc.createVariable('lon',np.float64,('lon',))
  time = data_nc.createVariable('time',np.float64,('time',))
  date = data_nc.createVariable('date',np.int32,('time',))
  datesec = data_nc.createVariable('datesec',np.int32,('time',))

  # Set variable units
  data_var.units = 'm/s'
  lat.units = 'degrees_north'
  lon.units = 'degrees_east'
  time.units = 'days since 1948-01-01 00:00:00'
  datesec.units = 's'

  # Set variable long names
  data_var.long_name = var
  lat.long_name = 'latitude'
  lon.long_name = 'longitude'
  time.long_name = 'observation time'
  date.long_name = 'current date as 8 digit integer (YYYYMMDD)'
  datesec.long_name = 'seconds to complete currents date'

  # Set other variable attributes
  time.calendar = 'noleap'
  time.note = 'time = 365 reset to 364.9999 so it is in same year as the preceeding data'
  data_var.level = '10m extracted'
  data_var.sampling_time = 'hourly'
  data_var.note = 'CFSR hourly data'

  # Set global attributes
  data_nc.title = 'CFSR hourly reanalysis time series'
  data_nc.Conventions = 'CF-1.0'
  data_nc.source = 'https://rda.ucar.edu/datasets/ds093.1/'
  data_nc.history = 'Interpolated from 0.5x0.5 degree CFSR product'

  # Set variables
  data_var[:,:,:] = data
  lat[:] = np.linspace(-90.0,90.0,nlat)
  lon[:] = np.linspace(-180.0,180.0,nlon)

  t = 0.0
  dt = 1.0/24.0
  for i in range(nsnaps):
    t = t + dt
    time[i] = t 

  t = 0.0
  dt = 3600
  for i in range(nsnaps):
    t = t + dt
    datesec[i] = t % 86400

  t = datetime.datetime(1948,1,1,0,0,0)
  for i in range(nsnaps): 
    t = t + datetime.timedelta(seconds=3600.0)
    date[i] = int(t.strftime('%Y%m%d'))
  
  
  data_nc.close()

##################################################################################################
##################################################################################################

def plot_interp_data(lon_data,lat_data,data,lon_grid,lat_grid,interp_data,var_label):

  levels = np.linspace(np.amin(data),np.amax(data),10)

  fig = plt.figure()
  ax0 = fig.add_subplot(2,1,1)
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

  parser = argparse.ArgumentParser()
  parser.add_argument('--plot',action='store_true')
  args = parser.parse_args()
  
  # Files to interpolate to/from
  grid_file = '/users/sbrus/scratch3/ACME/input_data/ocn/iaf/ncep.u_10.T62.1948.nc'
  data_file = '/users/sbrus/scratch4/WW3_testing/wind_data/june/wnd10mx0.5.gdas.200506.ww3.nc'

  # Interpolation of u and v velocities
  lon_grid,lat_grid,u_interp,lon_data,lat_data,u_data = interpolate_data_to_grid(grid_file,data_file,'U_GRD_L103')
  lon_grid,lat_grid,v_interp,lon_data,lat_data,v_data = interpolate_data_to_grid(grid_file,data_file,'V_GRD_L103')
 
  # Calculate and plot velocity magnitude
  if args.plot:
    for i in range(u_data.shape[0]):
      print i 
  
      data = np.sqrt(np.square(u_data[i,:,:]) + np.square(v_data[i,:,:]))
      interp_data = np.sqrt(np.square(u_interp[i,:,:]) + np.square(v_interp[i,:,:]))
      
      plot_interp_data(lon_data,lat_data,data,lon_grid,lat_grid,interp_data,'velocity magnitude')

  write_to_file('ncep.u_10.T62.1948.nc',u_interp,'u_10')
  write_to_file('ncep.v_10.T62.1948.nc',v_interp,'v_10')
