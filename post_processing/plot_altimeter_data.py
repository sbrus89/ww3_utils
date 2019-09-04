import netCDF4
import glob
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import datetime
import calendar
from skimage.measure import block_reduce
plt.switch_backend('agg')

write_nc = True

###############################################################################################
###############################################################################################

def read_field_ww3(f,variable='hs'):

  nc_file = netCDF4.Dataset(f,'r')

  ref_date = nc_file.variables['time'].getncattr('units').replace('days since ','')
  ref_date = datetime.datetime.strptime(ref_date,'%Y-%m-%d %H:%M:%S')
  output_time = nc_file.variables['time'][:]

  if len(output_time) > 1:
    print "Should be only one timestep per nc file. Check ww3_ounf.inp"
    raise SystemExit(0)

  date = ref_date + datetime.timedelta(days=output_time[0])
  output_date = date.strftime('%Y %m %d %H %M')

  lon  = nc_file.variables['longitude'][:]
  lat  = nc_file.variables['latitude'][:]

  var = nc_file.variables[variable][0,:,:]

  var = block_reduce(var, block_size=(2,2), func=np.mean)
  lon = lon[0::2]
  lat = lat[0::2]

  return lon,lat,var,output_date

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

def compute_altimeter_average(year_start,year_end,altimeter_files,ww3_files,month=None,season=None,year=None):

  swh_obs = np.zeros((180,360))
  nobs = np.zeros((180,360))

  swh_model = np.zeros((180,360))
  nmodel = 0.0

  arg_count = 0
  if month != None:
    arg_count = arg_count + 1
  #if season != None:
  #  arg_count = arg_count + 1
  if year != None:
    arg_count = arg_count + 1

  if arg_count == 0 or arg_count > 1 :
    raise SystemExit(0)

  if month:
    months = [month]
    year_adj = [0]
  #elif season == 'spring':
  #  months = [3,4,5]
  #  year_adj = [0,0,0]
  #elif season == 'summer':
  #  months = [6,7,8]
  #  year_adj = [0,0,0]
  #elif season == 'fall':
  #  months = [9,10,11]
  #  year_adj = [0,0,0]
  #elif season == 'winter':
  #  months = [12,1,2]
  #  year_adj = [-1,0,0]
  elif year:
    months   = [1,2,3,4,5,6,7,8,9,10,11,12]
    year_adj = [0,0,0,0,0,0,0,0,0, 0, 0, 0]
    
  for yr in range(year_start,year_end+1): 
    for i,mnth in enumerate(months):
      year = yr + year_adj[i]
  
      match = None
      for f in altimeter_files:
        if f.split('/')[-1].find(str(yr)) > 0:
          match = f
      if match:
        print match
        sig_wave_height,num_observations,lon,lat = read_altimeter_data(match)
      else:
        print('Error in finding altimeter files')
        raise SystemExit(0)

      # Accumulate SWH observations and counts
      swh_obs = swh_obs + sig_wave_height[mnth-1,0:-1,0:-1]
      nobs = nobs + num_observations[mnth-1,0:-1,0:-1]  
      lon_vec = lon[0:-1]
      lat_vec = lat[0:-1]

      month_start = datetime.datetime(year,mnth,1,0,0)
      month_end = datetime.datetime(year,mnth,calendar.monthrange(year,mnth)[1],23,59)

      for f in ww3_files:
        filename = f.split('/')[-1]

        file_date = filename.split('.')[1]
        try:
          file_datetime = datetime.datetime.strptime(file_date,'%Y%m%dT%HZ')
        except ValueError:
          continue

        if (file_datetime >= month_start) and (file_datetime <= month_end):
          print filename
          lon_vec,lat_vec,hs,output_date = read_field_ww3(f)
          swh_model = swh_model + hs
          nmodel = nmodel + 1.0

  swh_obs = np.divide(swh_obs,nobs.astype(np.float))
  swh_model = swh_model/nmodel

  return lon_vec,lat_vec,swh_obs,swh_model

###############################################################################
###############################################################################

def plot_altimeter_data(lon_vec,lat_vec,swh,filename,plot_type='field'):

  if plot_type == 'field':
    cmap = 'viridis'
    levels = np.linspace(0.0,10.0,50)
    label = 'Significant wave height (m)'
    ticks = []
  elif plot_type == 'difference':
    cmap = 'bwr'
    abs_max = np.nanmax(np.absolute(swh))
    levels = np.linspace(-abs_max,abs_max,100)
    ticks = [-abs_max, -0.5*abs_max, 0.0, 0.5*abs_max,abs_max]
    label = 'Percent difference in SWH'

  plt.figure(figsize=(7, 4.5))
  m = Basemap(projection='cyl',lon_0=0.0, llcrnrlat= -90,urcrnrlat=90,\
                               llcrnrlon=0,urcrnrlon=360,resolution='c')
  m.fillcontinents(color='tan',lake_color='lightblue')
  m.drawcoastlines()

  cf = plt.contourf(lon_vec,lat_vec,swh,levels,cmap=cmap)
  cb = plt.colorbar(cf,orientation='horizontal')
  cb.set_label(label)
  if ticks:
    cb.set_ticks(ticks)
  plt.tight_layout()
  plt.savefig(filename,bbox_inches='tight')
  plt.close()

###############################################################################
###############################################################################

def plot_comparison(lon_vec,lat_vec,swh_obs,swh_model,filename):

  # Create gridded latitude array
  Lon,Lat = np.meshgrid(lon_vec,lat_vec)
 
  # Get rid bad values
  idx = np.where(np.absolute(swh_model) > 1e10)
  swh_model[idx] = np.nan

  # Set max for x and y axes
  #swh_max = max(np.nanmax(swh_model),np.nanmax(swh_obs))
  swh_max = 9.0
 
  # Limit to only low latitudes  
  #idx = np.where(Lat > -40.0)
  #swh_model[idx] = np.nan

  fig = plt.figure()
  ax = fig.gca()
  sc = ax.scatter(swh_obs.ravel(),swh_model.ravel(),c=Lat.ravel(),vmin=-80.0,vmax=80.0)
  cb = plt.colorbar(sc,orientation='vertical')
  ax.plot([0.0,swh_max],[0.0,swh_max],'k')
  ax.set_xlim(0.0,swh_max)
  ax.set_ylim(0.0,swh_max)
  ax.set_xlabel('Observed')
  ax.set_ylabel('Modeled')
  cb.set_label('Degrees Latitude')
  fig.tight_layout()
  plt.savefig(filename,bbox_inches='tight')
  plt.close()
  

###############################################################################
###############################################################################

def create_plots(lon_vec,lat_vec,swh_obs,swh_model,altimeter_files,ww3_files,filename_id):

    

  if altimeter_files:
    filename = 'swh_obs_avg_'+str(year_start)+'-'+str(year_end)+filename_id+'.png'
    plot_altimeter_data(lon_vec,lat_vec,swh_obs,filename)
  if ww3_files:
    filename = 'swh_ww3_avg_'+str(year_start)+'-'+str(year_end)+filename_id+'.png'
    plot_altimeter_data(lon_vec,lat_vec,swh_model,filename)
  if altimeter_files and ww3_files:
    filename = 'swh_diff_avg_'+str(year_start)+'-'+str(year_end)+filename_id+'.png'
    diff = np.divide(swh_model-swh_obs,swh_obs)*100.0
    idx = np.where(np.absolute(diff) > 1e10)
    diff[idx] = np.nan
    plot_altimeter_data(lon_vec,lat_vec,diff,filename,plot_type='difference')
    filename = 'swh_comp_avg_'+str(year_start)+'-'+str(year_end)+filename_id+'.png'
    plot_comparison(lon_vec,lat_vec,swh_obs,swh_model,filename)

  if write_nc:
    filename = filename.split('.')[0]+'.nc'
    nc_file = netCDF4.Dataset(filename,'w')
    nc_file.createDimension('lon',lon_vec.size)
    nc_file.createDimension('lat',lat_vec.size)
    lon_nc = nc_file.createVariable('lon','f8',('lon',))
    lat_nc = nc_file.createVariable('lat','f8',('lat',))
    model_nc = nc_file.createVariable('swh_model','f8',('lat','lon'))
    obs_nc = nc_file.createVariable('swh_obs','f8',('lat','lon'))
    lon_nc[:] = lon_vec
    lat_nc[:] = lat_vec
    model_nc[:] = swh_model
    obs_nc[:] = swh_obs
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
    
  

###############################################################################
###############################################################################

if __name__ == '__main__':

    monthly = True 
    seasonal = False
    yearly = False 

    ww3_direcs = ['/users/sbrus/scratch4/WW3_CFSR_2000-2010/2002/run00/results/model_data/fields/',
                  '/users/sbrus/scratch4/WW3_CFSR_2000-2010/2003/run00/results/model_data/fields/']
    altimeter_direcs = ['/users/sbrus/scratch4/WW3_CFSR_2000-2010/2002/observed_data/altimeter/',
                        '/users/sbrus/scratch4/WW3_CFSR_2000-2010/2003/observed_data/altimeter/']

    # Find altimeter data files
    altimeter_files = []
    for direc in altimeter_direcs:
      altimeter_files.extend(glob.glob(direc+'*.nc'))
    altimeter_files.sort()
    print altimeter_files

    # Find model output files
    ww3_files = []
    for direc in ww3_direcs:
      ww3_files.extend(glob.glob(direc+'*.nc'))
    ww3_files.sort()

    # Get start and end years from filenames
    year_start = int(ww3_files[0].split('/')[-1].split('.')[1][0:4])
    year_end = int(ww3_files[-1].split('/')[-1].split('.')[1][0:4])-1

  
    if monthly == True:
      for mnth in range(1,13):
   
        lon_vec,lat_vec,swh_obs,swh_model = compute_altimeter_average(year_start,year_end,altimeter_files,ww3_files,month=mnth)
        create_plots(lon_vec,lat_vec,swh_obs,swh_model,altimeter_files,ww3_files,'_month'+str(mnth).zfill(2))        

    if seasonal == True:
      for season in ['winter','spring','summer','fall']:

        lon_vec,lat_vec,swh = compute_altimeter_average(year_start,year_end,altimeter_files,ww3_files,season=season)
        create_plots(lon_vec,lat_vec,swh_obs,swh_model,altimeter_files,ww3_files,'_'+season)        

    if yearly == True:
      for year in range(year_start,year_end+1):
 
        lon_vec,lat_vec,swh = compute_altimeter_average(year,year,altimeter_files,ww3_files,year=True)
        create_plots(lon_vec,lat_vec,swh_obs,swh_model,altimeter_files,ww3_files,'_year'+str(year).zfill(2))        
   
