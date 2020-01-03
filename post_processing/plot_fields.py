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
from scipy import interpolate
from mpl_toolkits.basemap import Basemap
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

  var = nc_file.variables[variable][0,:,:]
  nc_file.close()

  return lon,lat,var,output_date

###############################################################################################
###############################################################################################

def difference_fields(lon,lat,var,output_date,sign,filename,variable='hs'):

  lon2,lat2,var2,output_date2 = read_field_ww3(filename,variable)
  if var.shape == var2.shape:
    if sign == '-':
      var = var-var2
    else: 
      var = var+var2
  elif (var.shape[0] == var2.shape[0]/2) and (var.shape[1] == var2.shape[1]/2):
    if sign == '-':
      var = var-var2[::2,::2]
    else:
      var = var+var2[::2,::2]
  elif (var.shape[0] == var2.shape[0]*2) and (var.shape[1] == var2.shape[1]*2):
    if sign == '-':
      var = var[::2,::2]-var2
    else:
      var = var[::2,::2]+var2
    lon = lon2
    lat = lat2

  if output_date != output_date2:
    print('dates do not match')
    raise SystemExit(0)

  return lon,lat,var,output_date

###############################################################################################
###############################################################################################

def get_field_reanalysis(year,lon,lat,t,reanalysis_files,variable='ICEC_L1'):

  for f in reanalysis_files:
    filename = f.split('/')[-1]
    if filename.find(str(year)) > 0:

      nc_file = netCDF4.Dataset(f,'r')
    
      ref_date = nc_file.variables['time'].getncattr('units').replace('hours since ','').replace('.0 +0:00','')
      ref_date = datetime.datetime.strptime(ref_date,'%Y-%m-%d %H:%M:%S')
      data_time = nc_file.variables['time'][:]

      for i,time in enumerate(data_time):
        time_check = ref_date + datetime.timedelta(hours=int(time))

        if time_check == t:

          print("matching reanalysis data found")
 
          lon_data = nc_file.variables['lon'][:]
          lat_data = nc_file.variables['lat'][:]
          data = nc_file.variables[variable][i,:,:]
          nc_file.close()

          interp = interpolate.RegularGridInterpolator((lon_data,lat_data), data.T, bounds_error=False, fill_value=np.nan)  
           
          Lon,Lat = np.meshgrid(lon,lat)
          xypts = np.vstack((Lon.ravel(),Lat.ravel())).T

          var_interp = interp(xypts)
          var_interp = np.reshape(var_interp,Lat.shape)

          #fig = plt.figure(figsize=[18.0,9.0])
          #m = Basemap(projection='cyl',llcrnrlat= -90,urcrnrlat=90,\
          #                             llcrnrlon=0,urcrnrlon=360,resolution='c')
          #m.fillcontinents(color='tan',lake_color='lightblue')
          #m.drawcoastlines()
          #cf = plt.contourf(lon,lat,var_interp,100)
          #plt.title('interpolated ice')
          #cbar = plt.colorbar(cf,orientation='horizontal')
          #plt.savefig('interpolated_ice_'+time_check.strftime('%Y-%m-%d-%H')+'.png')
          #plt.close()

          return var_interp

###############################################################################################
###############################################################################################


def compute_metric(metric,year_start,year_end,solutions,month=None,year=None,season=None,reanalysis_files=None):

  arg_count = 0
  if month != None:
    arg_count = arg_count + 1
  if season != None:
    arg_count = arg_count + 1
  if year != None:
    arg_count = arg_count + 1

  if arg_count == 0 or arg_count > 1 :
    print('only one averaging period can be requested')
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
  elif year:
    months   = [1,2,3,4,5,6,7,8,9,10,11,12]
    year_adj = [0,0,0,0,0,0,0,0,0, 0, 0, 0]

  for i,run in enumerate(solutions):
    if i == 0:
      nfiles = len(run['files'])
    if nfiles != len(run['files']):
      print('number of files for solutions do not match')
      raise SystemExit(0)
 
  
  count = 0
  for year in range(year_start,year_end+1):
    for mnth in months:

      month_start = datetime.datetime(year,mnth,1,0,0)
      month_end = datetime.datetime(year,mnth,calendar.monthrange(year,mnth)[1],23,59)

      for i in range(nfiles):
         
        for j,run in enumerate(solutions):
          f = run['files'][i]

          filename = f.split('/')[-1]
                  
          file_date = filename.split('.')[1]
          file_datetime = datetime.datetime.strptime(file_date,'%Y%m%dT%HZ')
   
          in_range = False
          if (file_datetime >= month_start) and (file_datetime <= month_end):
            in_range = True
            print(run['sign'],f)
     
            # Read in and compute fields to plot
            if j == 0:
              lon,lat,var,output_date = read_field_ww3(f)
              if run['sign'] == '-':
                var = var*-1.0
            else:
              lon,lat,var,output_date = difference_fields(lon,lat,var,output_date,run['sign'],f)

            
        if in_range:

          if reanalysis_files:
            reanalysis_data = get_field_reanalysis(year,lon,lat,file_datetime,reanalysis_files,variable='ICEC_L1')
            reanalysis_data = ma.masked_outside(reanalysis_data,0.18,0.80) 

            var = ma.array(var,mask=reanalysis_data.mask)
            var = var.filled(fill_value=0.0)
        

          if metric == 'average':
            # Initialize average array
            if count == 0:
              count_array = np.zeros(var.shape)
              ones_array = np.ones(var.shape) 
              var_avg = np.zeros(var.shape)
            # Deal with mask from reanalysis
            if reanalysis_files:
              ones = ma.array(ones_array,mask=reanalysis_data.mask)
              ones = ones.filled(fill_value=0.0)
            else:
              ones = ones_array
            # Accumulate for average 
            count_array = count_array + ones 
            var_avg = var_avg + var

            #fig = plt.figure(figsize=[18.0,9.0])
            #m = Basemap(projection='cyl',llcrnrlat= -90,urcrnrlat=90,\
            #                             llcrnrlon=0,urcrnrlon=360,resolution='c')
            #m.fillcontinents(color='tan',lake_color='lightblue')
            #m.drawcoastlines()
            #cf = plt.contourf(lon,lat,var,100)
            #plt.title('ice mask')
            #cbar = plt.colorbar(cf,orientation='horizontal')
            #plt.savefig('ice_mask_'+file_datetime.strftime('%Y-%m-%d-%H')+'.png')
            #plt.close()

          elif metric == 'max':
            # Initialize max array
            if count == 0:
              var_avg = np.zeros_like(var)-1e36
            # Deal with masked array
            var = var.filled(fill_value=-1e36)
            var_avg = np.maximum(var_avg,var)
          elif metric == 'absmax':
            # Initialize max array
            if count == 0:
              var_avg = np.zeros_like(var)-1e36
            # Deal with masked array
            var = var.filled(fill_value=-1e36)
            var_avg = np.maximum(var_avg,np.absolute(var))
            print(np.amax(var_avg))
          elif metric == 'min':
            # Initialize min array
            if count == 0:
              var_avg = np.zeros_like(var)+1e36
            # Deal with masked array
            var = var.filled(fill_value=1e36)
            var_avg = np.minimum(var_avg,var)

          count = count + 1
          print("")
    
  # Compute average
  if metric == 'average':
    #var_avg = var_avg/float(count)
    var_avg = np.divide(var_avg,count_array)

  #return lon,lat,var_avg
  return lon,lat,var_avg

###############################################################################################
###############################################################################################

def plot_timesnaps(files,cfg):

  nruns,cmap,diff,run,symmetric_range = determine_plot_type(cfg)

  for f in files:
    filename = f.split('/')[-1]
    print(filename)
  
    # Read in and compute fields to plot
    if nruns > 1:
      lon,lat,var,output_date = difference_fields(f,cfg['model_direc'][1][1]+filename)
    else:
      lon,lat,var,output_date = read_field_ww3(f)

    # Plot timesnap
    if cfg["timesnaps"]:    
      title = 'Wave height'+diff+run+' on '+output_date
      cbar_label = 'wave height'+diff+' (m)'
      filename = run+'_'+'-'.join(output_date.split())+'.png'
      if "tsnap_range" in cfg:
        plot_field(lon,lat,var,cmap,title,cbar_label,filename,cfg['tsnap_range'][0],cfg['tsnap_range'][1])  
      else:
        plot_field(lon,lat,var,cmap,title,cbar_label,filename,symmetric_range=symmetric_range)  

###############################################################################################
###############################################################################################

def determine_plot_type(metric='average',run_list=[],nruns=0,range_min=None,range_max=None):
  
  # Find nc file names
  if nruns == 0:
    runs = [x[0] for x in run_list]
    nruns = len(runs)
  
  # Determine plot type based on number of runs
  if nruns > 2:
    print("Only 2 runs can be specified for difference plots")
    raise SystemExit(0)
  if nruns == 1:
    cmap = 'viridis'
    diff = ''
    symmetric_range = False
    run = run_list[0][0]
  elif nruns == 2 and not range_min and not range_max and metric == 'average': 
    cmap = 'bwr'
    diff = 'difference'
    symmetric_range = True
    run = ''
  elif nruns == 2 and metric == 'absmax':
    cmap = 'viridis'
    diff = 'difference'
    symmetric_range = False
    run = ''
  elif nruns == 2 and abs(range_min) == range_max and metric == 'average': 
    cmap = 'bwr'
    diff = 'difference'
    symmetric_range = True
    run = ''
  elif nruns == 2:
    cmap = 'viridis'
    diff = 'difference'
    symmetric_range = False
    run = ''

  return nruns,cmap,diff,run,symmetric_range

###############################################################################################
###############################################################################################

def plot_field(lon,lat,var,cmap,title,cbar_label,filename,range_min=None,range_max=None,symmetric_range=False,write_nc=False):

  idx = np.where(np.absolute(var) > 1e10)
  var[idx] = np.nan

  fig = plt.figure(figsize=[18.0,9.0])
  m = Basemap(projection='cyl',llcrnrlat= -90,urcrnrlat=-40,\
                               llcrnrlon=0,urcrnrlon=360,resolution='c')
  m.fillcontinents(color='tan',lake_color='lightblue')
  m.drawcoastlines()
  ticks = []
  if range_min != None and range_max != None:
    levels = np.linspace(range_min,range_max,100)
    if abs(range_min) == range_max:
      ticks = [-range_max, -0.5*range_max, 0.0, 0.5*range_max, range_max]
  elif symmetric_range:
    abs_max = np.nanmax(np.absolute(var))
    levels = np.linspace(-abs_max,abs_max,100)
    ticks = [-abs_max -0.5*abs_max, 0.0, 0.5*abs_max, abs_max]
    print(levels)
  else:
    levels = 100

   
  cf = plt.contourf(lon,lat,var,levels,cmap=cmap)
  plt.title(title)
  cbar = plt.colorbar(cf,orientation='horizontal')
  cbar.set_label(cbar_label)
  if ticks:
    cbar.set_ticks(ticks)
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
 
  monthly = True 
  yearly = False 
  overall = False 
 
  # Read in configuration file
  f = open(pwd+'/plot_fields.config')
  cfg = yaml.load(f,Loader=yaml.Loader)
  pprint.pprint(cfg)
 
  # Set colorbar range if specified
  avg_range_min = None 
  avg_range_max = None
  if 'avg_range' in cfg:
    avg_range_min = cfg['avg_range'][0]
    avg_range_max = cfg['avg_range'][1]

  if 'metric' in cfg:
    metric = cfg['metric']
  else:
    metric = 'average'

  if 'reanalysis_files' in cfg:
    reanalysis_files = cfg['reanalysis_files']
  else:
    reanalysis_files = []

  # Determine if single solution or difference of two is requested
  nitems = len(cfg['model_runs'])
  if nitems > 1:
    nruns = 2
  else:
    nruns = 1

  nruns,cmap,diff,run,symmetric_range = determine_plot_type(metric=metric,nruns=nruns,range_min=avg_range_min,range_max=avg_range_max)

  solutions = [] 
  for j,item in enumerate(cfg['model_runs']):   
    
    # Get field output files
    run_name = item[0] 
    solutions.append({'name':run_name,'sign':item[1],'files':[]})
    for output_direc in item[2:]:
      solutions[j]['files'].extend(glob.glob(output_direc+'*.nc'))
    solutions[j]['files'].sort()
  
    # Get start and end years from filenames
    start = int(solutions[j]['files'][0].split('/')[-1].split('.')[1][0:4])
    end = int(solutions[j]['files'][-1].split('/')[-1].split('.')[1][0:4])-1
    if j == 0:
      start_year = start
      end_year = end
      print(start_year)
      print(end_year)
    else:
      if (start != start_year) or (end != end_year):
        raise SystemExit(0)

  pprint.pprint(solutions)

  # Set averaging period range
  if monthly:
    avg_periods = range(1,13)
    yearly = None
  elif yearly:
    avg_periods = range(start_year,end_year)
    mnth = None
  elif overall:
    avg_periods = range(1)
    yearly = True
    mnth = None
  
      
  for it,i in enumerate(avg_periods): 
    print(i)

    if monthly:
      mnth = i
      year_start = start_year
      year_end = end_year
    elif overall:
      mnth = None
      year_start = start_year
      year_end = end_year
    elif yearly:
      mnth = None
      year_start = i
      year_end = i

    if monthly:
      saved_results = metric+'_'+diff+run+'_month'+str(mnth).zfill(2)+'_year'+str(year_start)+'-'+str(year_end)+'.nc'
    elif yearly:
      saved_results = metric+'_'+diff+run+'_year'+str(year_start)+'-'+str(year_end)+'.nc'

    write_nc = True
    if os.path.exists(saved_results):
      nc_file = netCDF4.Dataset(saved_results,'r')    
      lon = nc_file.variables['lon'][:]
      lat = nc_file.variables['lat'][:]
      var = nc_file.variables['var'][:,:]
      nc_file.close()
      write_nc = False
    else:
      lon,lat,var = compute_metric(metric,year_start,year_end,solutions,month=mnth,year=yearly,reanalysis_files=reanalysis_files)

    if it == 0:
      print(var.shape)
      lat_trend = np.zeros((lat.size,len(avg_periods)))
    lat_trend[:,it] = np.mean(var,axis=1)
      

    # Plot average
    cbar_label = 'wave height '+diff+' (m)'
    if monthly:
      title = metric +' wave height '+diff+run+' for month '+str(mnth).zfill(2)+' years '+str(year_start)+'-'+str(year_end)
      filename = metric+'_'+diff+run+'_month'+str(mnth).zfill(2)+'_year'+str(year_start)+'-'+str(year_end)+'.png'
    elif yearly:
      title = metric+' wave height '+diff+run+' years '+str(year_start)+'-'+str(year_end)
      filename = metric+'_'+diff+run+'_year'+str(year_start)+'-'+str(year_end)+'.png'
    plot_field(lon,lat,var,cmap,title,cbar_label,filename,avg_range_min,avg_range_max,symmetric_range,write_nc)  


  title = "average variation"
  if monthly:
    filename = "monthly_variation.png"
  if yearly:
    filename = "yearly_variation.png"
  fig = plt.figure(figsize=[18.0,9.0])
  levels = 100
  cf = plt.contourf(avg_periods,lat,lat_trend,levels,cmap=cmap)
  plt.title(title)
  cbar = plt.colorbar(cf,orientation='horizontal')
  cbar.set_label(cbar_label)
  plt.savefig(filename)
  plt.close() 

  # Move plots to their own directory
  if not os.path.exists(cfg["plot_direc"]):
    subprocess.call(['mkdir','-p',cfg["plot_direc"]])
  subprocess.call('mv *.png '+cfg["plot_direc"],shell=True)    


