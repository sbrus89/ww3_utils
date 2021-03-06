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
    print 'Interpolating '+var+': '+str(i)

    # Get data to interpolate
    data[i,:,0:-1] = np.flipud(data_nc.variables[var][i,:,:])
    data[i,:,-1] = data[i,:,0]

    # Interpolate data onto new grid
    interpolator = interpolate.RegularGridInterpolator((lon_data,lat_data),data[i,:,:].T,bounds_error=False,fill_value=0.0)
    interp_data[i,:] = interpolator(grid_points)

  # Deal with time
  ref_date = data_nc.variables['time'].getncattr('units').replace('hours since ','').replace('.0 +0:00','')
  ref_date = datetime.datetime.strptime(ref_date,'%Y-%m-%d %H:%M:%S')
  xtime = []
  for t in time:
    date = ref_date + datetime.timedelta(hours=np.float64(t))
    xtime.append(date.strftime('%Y-%m-%d_%H:%M:%S'+45*' '))
  xtime = np.asarray(xtime)

  return lon_grid,lat_grid,interp_data,lon_data,lat_data,data,xtime

##################################################################################################
##################################################################################################

def write_to_file(filename,data,var,xtime):

  if os.path.isfile(filename):
    data_nc = netCDF4.Dataset(filename,'a', format='NETCDF3_64BIT_OFFSET')
  else:
    data_nc = netCDF4.Dataset(filename,'w', format='NETCDF3_64BIT_OFFSET')

    # Find dimesions
    ncells = data.shape[1]
    nsnaps = data.shape[0]

    # Declare dimensions
    data_nc.createDimension('nCells',ncells)
    data_nc.createDimension('StrLen',64)
    data_nc.createDimension('Time',None)

    # Create time variable    
    time = data_nc.createVariable('xtime','S1',('Time','StrLen')) 
    time[:] = netCDF4.stringtochar(xtime) 

  # Declear variables
  time = data_nc.dimensions['Time'].name
  ncells = data_nc.dimensions['nCells'].name
  data_var = data_nc.createVariable(var,np.float64,(time,ncells))

  # Set variables
  data_var[:,:] = data[:,:]
  data_nc.close()

##################################################################################################
##################################################################################################

def plot_interp_data(lon_data,lat_data,data,lon_grid,lat_grid,interp_data,var_label,var_abrev,time):

  levels = np.linspace(np.amin(data),np.amax(data),10)

  # Plot data
  fig = plt.figure()
  ax0 = fig.add_subplot(2,1,1)
  cf = ax0.contourf(lon_data,lat_data,data,levels=levels)
  m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
                               llcrnrlon=0,urcrnrlon=360,resolution='c')
  m.fillcontinents(color='tan',lake_color='lightblue')
  m.drawcoastlines()
  ax0.set_title('data '+time.strip())
  cbar = fig.colorbar(cf,ax=ax0)
  cbar.set_label(var_label)

  # Plot interpolated data
  ax1 = fig.add_subplot(2,1,2)
  cf = ax1.tricontourf(lon_grid,lat_grid,interp_data,levels=levels)
  m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
                               llcrnrlon=0,urcrnrlon=360,resolution='c')
  m.fillcontinents(color='tan',lake_color='lightblue')
  m.drawcoastlines()
  ax1.set_title('interpolated data '+time.strip())
  cbar = fig.colorbar(cf,ax=ax1)
  cbar.set_label(var_label)

  # Save figure
  fig.tight_layout()
  fig.savefig(var_abrev+'_'+str(i).zfill(4)+'.png',box_inches='tight')
  plt.close()

##################################################################################################
##################################################################################################

if __name__ == '__main__':

  parser = argparse.ArgumentParser()
  parser.add_argument('--plot',action='store_true')
  args = parser.parse_args()
 
  nplot = 10
 
  # Files to interpolate to/from
  grid_file = '/users/sbrus/scratch4/MPAS-O_testing/ocean/hurricane/USDEQU120cr10rr1/build_mesh/culled_mesh/culled_mesh.nc'
  wind_file = '/users/sbrus/scratch4/MPAS-O_testing/ocean/hurricane/USDEQU120cr10rr1/time_varying_forcing/wind_data/wnd10m.nc'
  pres_file = '/users/sbrus/scratch4/MPAS-O_testing/ocean/hurricane/USDEQU120cr10rr1/time_varying_forcing/wind_data/prmsl.nc'
  forcing_file = 'atmospheric_forcing.nc'

  # Interpolation of u and v velocities
  lon_grid,lat_grid,u_interp,lon_data,lat_data,u_data,xtime = interpolate_data_to_grid(grid_file,wind_file,'U_GRD_L103')
  lon_grid,lat_grid,v_interp,lon_data,lat_data,v_data,xtime = interpolate_data_to_grid(grid_file,wind_file,'V_GRD_L103')
 
  # Calculate and plot velocity magnitude
  if args.plot:
    for i in range(u_data.shape[0]):
      if i % nplot == 0:

        print 'Plotting vel: '+str(i)
  
        data = np.sqrt(np.square(u_data[i,:,:]) + np.square(v_data[i,:,:]))
        interp_data = np.sqrt(np.square(u_interp[i,:]) + np.square(v_interp[i,:]))
      
        plot_interp_data(lon_data,lat_data,data,lon_grid,lat_grid,interp_data,'velocity magnitude','vel',xtime[i])

  # Interpolation of atmospheric pressure
  lon_grid,lat_grid,p_interp,lon_data,lat_data,p_data,xtime = interpolate_data_to_grid(grid_file,pres_file,'PRMSL_L101')
  
  # Plot atmopheric pressure
  if args.plot:
    for i in range(p_data.shape[0]):
      if i % nplot == 0:

        print 'Plotting pres: '+str(i)

        plot_interp_data(lon_data,lat_data,p_data[i,:,:],lon_grid,lat_grid,p_interp[i,:],'atmospheric pressure','pres',xtime[i])
  
  # Write to NetCDF file
  subprocess.call(['rm',forcing_file])
  write_to_file(forcing_file,u_interp,'windSpeedU',xtime)
  write_to_file(forcing_file,v_interp,'windSpeedV',xtime)
  write_to_file(forcing_file,p_interp,'atmosPressure',xtime)

