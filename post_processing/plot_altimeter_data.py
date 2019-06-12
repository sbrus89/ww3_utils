import netCDF4
import glob
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
import numpy as np
import datetime

###############################################################################
###############################################################################

def read_altimeter_data(filename):

  nc_file = netCDF4.Dataset(filename,'r')
  swh = nc_file.variables['swh_sum'][:]
  nobs = nc_file.variables['nobs'][:]
  lon = nc_file.variables['lon'][:]
  lat = nc_file.variables['lat'][:]
  nc_file.close()

  return swh,nobs,lon,lat

###############################################################################
###############################################################################

def compute_altimeter_average(year_start,year_end,month=None,season=None):

  swh = np.zeros((181,361))
  nobs = np.zeros((181,361))

  if month != None and season != None:
    raise SystemExit(0)

  if month:
    months = [month]
    year_adj = [0]
  elif season == 'spring':
    months = [3,4,5]
    year_adj = [0,0,0]
  elif season == 'summer':
    months = [6,7,8]
    year_adj = [0,0,0]
  elif season == 'fall':
    months = [9,10,11]
    year_adj = [0,0,0]
  elif season == 'winter':
    months = [12,1,2]
    year_adj = [-1,0,0]
    
  for yr in range(year_start,year_end+1): 
    for i,mnth in enumerate(months):
      year = yr + year_adj[i]

      filename = 'accumulated_monthly_'+str(yr)+'.nc'
      sig_wave_height,num_observations,lon_vec,lat_vec = read_altimeter_data(filename)

      # Accumulate SWH observations and counts
      swh = swh + sig_wave_height[mnth-1,:,:]
      nobs = nobs + num_observations[mnth-1,:,:]  


  swh = np.divide(swh,nobs.astype(np.float))

  return lon_vec,lat_vec,swh

###############################################################################
###############################################################################

def plot_altimeter_data(lon_vec,lat_vec,swh,filename):

  plt.figure(figsize=(14, 6))
  m = Basemap(projection='cyl',lon_0=0.0, llcrnrlat= -90,urcrnrlat=90,\
                               llcrnrlon=0,urcrnrlon=360,resolution='c')
  m.fillcontinents(color='tan',lake_color='lightblue')
  m.drawcoastlines()

  levels = np.linspace(0.0,10.0,50)
  cf = plt.contourf(lon_vec,lat_vec,swh,levels)
  plt.colorbar(cf,orientation='horizontal')
  plt.savefig(filename)

###############################################################################
###############################################################################

if __name__ == '__main__':

    monthly = True 
    seasonal = True

    year_start = 1999
    year_end   = 2010 
  
    if monthly == True:
      for mnth in range(1,13):
   
        lon_vec,lat_vec,swh = compute_altimeter_average(year_start,year_end,month=mnth)
        filename = 'swh_avg_'+str(year_start)+'-'+str(year_end)+'_month'+str(mnth).zfill(2)+'.png'
        plot_altimeter_data(lon_vec,lat_vec,swh,filename)

    if seasonal == True:
      for season in ['winter','spring','summer','fall']:

        lon_vec,lat_vec,swh = compute_altimeter_average(year_start,year_end,season=season)
        filename = 'swh_avg_'+str(year_start)+'-'+str(year_end)+'_'+season+'.png'
        plot_altimeter_data(lon_vec,lat_vec,swh,filename)

