import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import netCDF4 
import glob
import datetime
import calendar
import os
import yaml
import pprint
import subprocess
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy import interpolate
plt.switch_backend('agg')


###############################################################################################
###############################################################################################

def read_field_ww3(f,variable='hs'):

  nc_file = netCDF4.Dataset(f,'r')

  ref_date = nc_file.variables['time'].getncattr('units').replace('days since ','')
  ref_date = datetime.datetime.strptime(ref_date,'%Y-%m-%d %H:%M:%S')
  output_time = nc_file.variables['time'][:]

  if len(output_time) > 1:
    print("Should be only one timestep per nc file. Check ww3_ounf.inp")
    raise SystemExit(0)

  date = ref_date + datetime.timedelta(days=output_time[0])
  output_date = date.strftime('%Y %m %d %H %M')

  lon  = nc_file.variables['longitude'][:] 
  lat  = nc_file.variables['latitude'][:]

  var = nc_file.variables[variable][0,:]
  nc_file.close()

  return lon,lat,var,output_date

###############################################################################################
###############################################################################################

def plot_timesnaps(files,cfg):

  for f in files:
    filename = f.split('/')[-1]
    print(filename)
  
    # Read in and compute fields to plot
    lon,lat,var,output_date = read_field_ww3(f)

    # Plot timesnap
    cmap = 'viridis'
    title = 'Wave height on '+output_date
    cbar_label = 'wave height (m)'
    filename = '-'.join(output_date.split())+'.png'
    plot_field(lon,lat,var,cmap,title,cbar_label,filename)  

###############################################################################################
###############################################################################################

def plot_field(lon,lat,var,cmap,title,cbar_label,filename,write_nc=False):

  idx = np.where(np.absolute(var) > 1e10)
  var[idx] = np.nan

  fig = plt.figure(figsize=[18.0,9.0])
  ax = fig.add_subplot(1,1,1, projection=ccrs.PlateCarree())
  levels = 100
   
  cf = plt.tricontourf(lon,lat,var,levels,cmap=cmap,  transform=ccrs.PlateCarree())
  ax.set_extent([-180.0,180.0,-90.0,90.0], crs=ccrs.PlateCarree())
  ax.add_feature(cfeature.LAND, zorder=100)
  ax.add_feature(cfeature.LAKES, alpha=0.5, zorder=101)
  ax.add_feature(cfeature.COASTLINE, zorder=101)
  ax.set_title(title)
  cbar = plt.colorbar(cf,orientation='horizontal')
  cbar.set_label(cbar_label)
  plt.savefig(filename)
  plt.close() 

  if write_nc:
    filename = filename.split('.')[0]+'.nc'
    nc_file = netCDF4.Dataset(filename,'w')
    nc_file.createDimension('lon',lon.size)
    nc_file.createDimension('lat',lat.size)
    lon_nc = nc_file.createVariable('lon','f8',('lon',))
    lat_nc = nc_file.createVariable('lat','f8',('lat',))
    var_nc = nc_file.createVariable('var','f8',('lat','lon'))
    lon_nc[:] = lon
    lat_nc[:] = lat
    var_nc[:] = var
    lon_nc.units = "degree_east"
    lon_nc.long_name = "longitude"
    lon_nc.standard_name = "longitude"
    lon_nc.valid_min = -180.0
    lon_nc.valid_max = 360.0
    lon_nc.axis = "N"
    lat_nc.units = "degree_north"
    lat_nc.long_name = "latitude"
    lat_nc.standard_name = "latitude"
    lat_nc.valid_min = -90.0
    lat_nc.valid_max = 180.0
    lat_nc.axis = "Y"
    nc_file.close()


###############################################################################################
###############################################################################################

if __name__ == '__main__':

  pwd = os.getcwd()
 
  # Read in configuration file
  f = open(pwd+'/plot_fields.config')
  cfg = yaml.load(f,Loader=yaml.Loader)
  pprint.pprint(cfg)

  files = glob.glob(cfg['model_direc']+'*.nc')

  plot_timesnaps(files,cfg)

  # Move plots to their own directory
  if not os.path.exists(cfg["plot_direc"]):
    subprocess.call(['mkdir','-p',cfg["plot_direc"]])
  subprocess.call('mv *.png '+cfg["plot_direc"],shell=True)    


