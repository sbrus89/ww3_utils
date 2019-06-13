import netCDF4 
import matplotlib.pyplot as plt

from mpl_toolkits.basemap import Basemap
import numpy as np
import glob
import datetime
import calendar
import os
import yaml
import pprint
import subprocess
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
    print "Should be only one timestep per nc file. Check ww3_ounf.inp"
    raise SystemExit(0)

  date = ref_date + datetime.timedelta(days=output_time[0])
  output_date = date.strftime('%Y %m %d %H %M')

  lon  = nc_file.variables['longitude'][:] 
  lat  = nc_file.variables['latitude'][:]

  var = nc_file.variables[variable][0,:,:]

  return lon,lat,var,output_date

###############################################################################################
###############################################################################################

def difference_fields(file1,file2,variable='hs'):

  lon,lat,var,output_date = read_field_ww3(file1,variable)
  lon2,lat2,var2,output_date2 = read_field_ww3(file2,variable)
  if var.shape == var2.shape:
    var = var-var2
  elif (var.shape[0] == var2.shape[0]/2) and (var.shape[1] == var2.shape[1]/2):
    var = var-var2[::2,::2]
  elif (var.shape[0] == var2.shape[0]*2) and (var.shape[1] == var2.shape[1]*2):
    var = var[::2,::2]-var2
    lon = lon2
    lat = lat2
  if output_date != output_date2:
    print 'dates do not match'
    raise SystemExit(0)

  return lon,lat,var,output_date

###############################################################################################
###############################################################################################

def compute_average(year_start,year_end,files,files2=None,month=None,year=None,season=None):

  arg_count = 0
  if month != None:
    arg_count = arg_count + 1
  if season != None:
    arg_count = arg_count + 1
  if year != None:
    arg_count = arg_count + 1

  if arg_count == 0 or arg_count > 1 :
    print 'only one averaging period can be requested'
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
 
  
  count = 0
  for year in range(year_start,year_end+1):
    for mnth in months:

      month_start = datetime.datetime(year,mnth,1,0,0)
      month_end = datetime.datetime(year,mnth,calendar.monthrange(year,mnth)[1],23,59)

      for i,f in enumerate(files):
        filename = f.split('/')[-1]
                
        file_date = filename.split('.')[1]
        file_datetime = datetime.datetime.strptime(file_date,'%Y%m%dT%HZ')
   
        if (file_datetime >= month_start) and (file_datetime <= month_end):
          print filename
     
          # Read in and compute fields to plot
          if files2:
            lon,lat,var,output_date = difference_fields(f,files2[i])
          else:
            lon,lat,var,output_date = read_field_ww3(f)
    
          # Initialize average array
          if count == 0:
            var_avg = np.zeros(var.shape)

          # Accumulate for average
          var_avg = var_avg + var
          count = count + 1
    
  # Compute average
  var_avg = var_avg/float(count)

  return lon,lat,var_avg

###############################################################################################
###############################################################################################

def plot_timesnaps(files,cfg):

  nruns,cmap,diff,run,symmetric_range = determine_plot_type(cfg)

  for f in files:
    filename = f.split('/')[-1]
    print filename
  
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

def determine_plot_type(cfg):
  
  # Find nc file names
  runs = [x[0] for x in cfg['model_direc']]
  nruns = len(runs)
  
  # Determine plot type based on number of runs
  if nruns > 2:
    print "Only 2 runs can be specified for difference plots"
    raise SystemExit(0)
  if nruns == 1:
    cmap = 'viridis'
    diff = ''
    symmetric_range = False
    run = cfg['model_direc'][0][0]
  elif nruns == 2:
    cmap = 'bwr'
    diff = 'difference'
    symmetric_range = True
    run = ''

  return nruns,cmap,diff,run,symmetric_range

###############################################################################################
###############################################################################################

def plot_field(lon,lat,var,cmap,title,cbar_label,filename,range_min=None,range_max=None,symmetric_range=False):

  fig = plt.figure()
  m = Basemap(projection='cyl',llcrnrlat= -90,urcrnrlat=90,\
                               llcrnrlon=0,urcrnrlon=360,resolution='c')
  m.fillcontinents(color='tan',lake_color='lightblue')
  m.drawcoastlines()
  if range_min and range_max:
    levels = np.linspace(range_min,range_max,100)
  elif symmetric_range:
    abs_max = np.amax(np.absolute(var))
    levels = np.linspace(-abs_max,abs_max,100)
  else:
    levels = 100
   
  cf = plt.contourf(lon,lat,var,levels,cmap=cmap)
  plt.title(title)
  cbar = plt.colorbar(cf,orientation='horizontal')
  cbar.set_label(cbar_label)
  plt.savefig(filename)
  plt.close() 

###############################################################################################
###############################################################################################

if __name__ == '__main__':

  pwd = os.getcwd()
 
  monthly = False 
  yearly = False 
  overall = True
 
  # Read in configuration file
  f = open(pwd+'/plot_fields.config')
  cfg = yaml.load(f)
  pprint.pprint(cfg)
  nruns,cmap,diff,run,symmetric_range = determine_plot_type(cfg)
  
  # Get field output files
  files = []
  for output_direc in cfg['model_direc'][0][1:]:
    files.extend(glob.glob(output_direc+'*.nc'))
  files.sort()

  files2 = None
  if nruns > 1:
    files2 = []
    for output_direc in cfg['model_direc'][1][1:]:
      files2.extend(glob.glob(output_direc+'*.nc'))
    files2.sort()

  # Get start and end years from filenames
  start_year = int(files[0].split('/')[-1].split('.')[1][0:4])
  end_year = int(files[-1].split('/')[-1].split('.')[1][0:4])-1
  print start_year
  print end_year


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

  avg_range_min = None
  avg_range_max = None 
  if 'avg_range' in cfg:
    avg_range_min = cfg['avg_range'][0]
    avg_range_max = cfg['avg_range'][1]
  
    

  
  for i in avg_periods: 
    print i

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
     

    lon,lat,var_avg = compute_average(year_start,year_end,files,files2,month=mnth,year=yearly)
 
    # Plot average
    cbar_label = 'wave height '+diff+' (m)'
    if monthly:
      title = 'Average wave height '+diff+run+' for month '+str(mnth).zfill(2)+' years '+str(year_start)+'-'+str(year_end)
      filename = 'average_'+diff+run+'_month'+str(mnth).zfill(2)+'_year'+str(year_start)+'-'+str(year_end)+'.png'
    elif yearly:
      title = 'Average wave height '+diff+run+' years '+str(year_start)+'-'+str(year_end)
      filename = 'average_'+diff+run+'_year'+str(year_start)+'-'+str(year_end)+'.png'
    plot_field(lon,lat,var_avg,cmap,title,cbar_label,filename,avg_range_min,avg_range_max,symmetric_range)  


  # Move plots to their own directory
  if not os.path.exists(cfg["plot_direc"]):
    subprocess.call(['mkdir','-p',cfg["plot_direc"]])
  subprocess.call('mv *.png '+cfg["plot_direc"],shell=True)    


