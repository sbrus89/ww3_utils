import netCDF4 
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.gridspec as gridspec
import numpy as np
import glob
import pprint
import datetime
import os
import yaml
import pprint
import subprocess
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy import interpolate

plt.switch_backend('agg')
np.seterr(divide='ignore', invalid='ignore')
cartopy.config['pre_existing_data_dir'] = \
        os.getenv('CARTOPY_DIR', cartopy.config.get('pre_existing_data_dir'))

#--------------------------
# Define variables to plot
#--------------------------
variables = {
            'hs'  :{'obs_col' : 8,
                     'fill_val': 99.00,
                     'recip'   : False,
                     'title'   : 'Significant wave height',
                     'label'   : 'wave height',
                     'units'   : 'm',
                     'aka'     : ['HS']},
             'th1p':{'obs_col' : 11, 
                     'fill_val': 999,
                     'recip'   : False,
                     'title'   : 'Dominant wave direction',
                     'label'   : 'wave direction',
                     'units'   : 'deg',
                     'aka'     : ['dp']},
             'fp'  :{'obs_col' : 9,
                     'fill_val': 99.0,
                     'recip'   : True,
                     'title'   : 'Dominant wave period',
                     'label'   : 'wave period',
                     'units'   : 's',
                     'aka'     : ['FP']},
             'wnd' :{'obs_col' : 6,
                     'fill_val': 99.0,
                     'recip'   : False,
                     'title'   : 'Wind speed',
                     'label'   : 'wind speed',
                     'units'   : 'm/s',
                     'aka'     : ['uwnd','vwnd']},
             'wnddir':{'obs_col' : 5,
                       'fill_val': 999,
                       'recip'   : False,
                       'title'   : 'Wind direction',
                       'label'   : 'wind direction',
                       'units'   : 'deg',
                       'aka'     : ['uwnd','vwnd']},
             #'ssh  ':{'obs_col' : 5,
             #          'fill_val': 99.0,
             #          'recip'   : False,
             #          'title'   : 'Water Level',
             #          'units'   : 'm',
             #          'aka'     : ['SSH']}
              }

################################################################################################
################################################################################################

def read_point_files(data_files,variables,start_date):

  for j,files in enumerate(data_files):
    for i,name in enumerate(files):
      print(name.split('/')[-1])
    
      nc_file = netCDF4.Dataset(name,'r')
    
      # Initializations for first iteration of run
      if j == 0 and i == 0:
          # Station information
          stations = {}
          stations['name'] = netCDF4.chartostring(nc_file.variables['station_name'][:,:]).tolist()
          stations['lon'] = np.squeeze(nc_file.variables['longitude'][:])
          stations['lat'] = np.squeeze(nc_file.variables['latitude'][:])
          nstations = len(stations['name'])
  
          # Model output data (nfiles = number of timesnaps, 1 timesnap per file)
          nfiles = len(files)
          data = {}
          for var in variables:
            data[var] = np.empty((nfiles,nstations))
    
          # Time information
          data['time'] = np.empty(nfiles)
          data['ref_date'] = nc_file.variables['time'].getncattr('units').replace('days since ','')
          data['ref_date'] = datetime.datetime.strptime(data['ref_date'],'%Y-%m-%d %H:%M:%S')
      
      # Get time and output variables
      if j == 0:
        data['time'][i] = nc_file.variables['time'][:]
      for var in variables:
        if var in nc_file.variables:
          print(var)
          data[var][i][:] = nc_file.variables[var][:]

  data['datetime'], data['date'] = output_time_to_date(data['time'],data['ref_date'],start_date)

  return data, stations

################################################################################################
################################################################################################

def interpolate_stations_from_fields(data_files,variables,station_file,start_date):

  stations = read_station_file(station_file)
  nstations = len(stations['name'])
  sta_pts = np.column_stack((stations['lon'],stations['lat']))

  t = 0.0
  for j,files in enumerate(data_files):
    for i,name in enumerate(files):
      print( name.split('/')[-1])

      nc_file = netCDF4.Dataset(name,'r')

      # Initializations
      if j == 0 and i == 0:
        
          # Model output data
          nfiles = len(data_files[0])
          data = {}
          for var in variables:
            data[var] = np.zeros((nfiles,nstations))
            #data[var].fill(np.nan)
  
          # Time information
          data['time'] = np.empty(nfiles)
          if 'time' in nc_file.variables:
            data['ref_date'] = nc_file.variables['time'].getncattr('units').replace('days since ','')
            data['ref_date'] = datetime.datetime.strptime(data['ref_date'],'%Y-%m-%d %H:%M:%S')
          else:
            # (for E3SM WW3 output)
            data['ref_date'] = datetime.datetime.strptime(cfg['ref_date'], '%Y-%m-%d %H:%M:%S')

          # Determine if unstructured
          unstructured = False
          if 'node' in nc_file.dimensions:
            unstructured = True
        
      # Get field grid points
      if 'longitude' in nc_file.variables:
        lon = nc_file.variables['longitude'][:]
        lat = nc_file.variables['latitude'][:]
      else:
        # (for E3SM WW3 output)
        lon  = np.linspace(cfg['lon_range'][0],cfg['lon_range'][1],nc_file.dimensions['NX'].size)
        lat  = np.linspace(cfg['lat_range'][0],cfg['lat_range'][1],nc_file.dimensions['NY'].size)
      idx = np.where(lon > 180.0)
      lon[idx] = lon[idx] - 360.0

      if unstructured == False:
        Lon,Lat = np.meshgrid(lon,lat)
      else:
        Lon = lon
        Lat = lat
        

      # Get time
      if j == 0:
        if 'time' in nc_file.variables:
          data['time'][i] = nc_file.variables['time'][:]
        else:
          # (for E3SM WW3 output)
          data['time'][i] = t
          t = t + cfg["out_int"] 

      # Get variables
      for var in variables:
        if var in nc_file.variables or any(x in nc_file.variables for x in variables[var]['aka']):
          print( var )     

          # Account for alernate variable names
          if var not in nc_file.variables:
            ls = [v for v in variables[var]['aka'] if v in nc_file.variables]
            if len(ls) == 1:                     
              ncvar = ls[0]                       # scalar fields 
              field = read_field(nc_file,ncvar)   # with alternate names
            if len(ls) == 2:
              ls.sort()                           # vector
              ncvar1 = ls[0]                      # fields
              field1 = read_field(nc_file,ncvar1)
              ncvar2 = ls[1]
              field2 = read_field(nc_file,ncvar2)
              if variables[var]['units'] == 'm/s':                       # magnitude
                field = np.sqrt(np.square(field1)+np.square(field2))  
              elif variables[var]['units'] == 'deg':                     # direction
                field = np.arctan2(field2,field1)
                field = np.degrees(field)
                #field = -field+360.0
          else:
            field = read_field(nc_file,var)      # scalar fields

          # Interpolate station values from field
          x = Lon[~field.mask]
          y = Lat[~field.mask]
          field = field[~field.mask]
          f = field.ravel()
          x = x.ravel()
          y = y.ravel()
          grid_pts = np.column_stack((x,y))
          if grid_pts.size > 0:
            interp = interpolate.LinearNDInterpolator(grid_pts,f)
            data[var][i][:] = interp(sta_pts)[:]

  data['datetime'], data['date'] = output_time_to_date(data['time'],data['ref_date'],start_date)

  return data, stations
################################################################################################
################################################################################################

def read_field(nc_file,var):

  if 'time' in nc_file.variables:
    if 'node' in nc_file.variables[var].dimensions:
      field = nc_file.variables[var][0,:]
    else:
      field = nc_file.variables[var][0,:,:]
  else:
    field = nc_file.variables[var][:,:]

  return field

################################################################################################
################################################################################################

def read_station_file(station_file):

  # Read in stations names and location
  f = open(station_file)
  lines = f.read().splitlines()
  stations = {}
  stations['name'] = []
  stations['lon'] = []
  stations['lat'] = []
  for sta in lines:
    val = sta.split()
    stations['name'].append(val[2].strip("'"))
    stations['lon'].append(float(val[0]))
    stations['lat'].append(float(val[1]))
  nstations = len(stations['name'])
  stations['lon'] = np.asarray(stations['lon'])
  stations['lat'] = np.asarray(stations['lat'])

  return stations 

################################################################################################
################################################################################################

def read_station_data(obs_file,min_date,max_date,variables,output_datetime=None):

  # Initialize variable for observation data
  obs_data = {}
  for var in variables:
    obs_data[var] = []  
  obs_data['datetime'] = []
 
  # Open the observation file
  f = open(obs_file)
  obs = f.read().splitlines()

  # Determine file format
  frmt = '%Y %m %d %H %M'
  if obs[0].find('YYYY MM DD hh mm') >= 0:
    obs_frmt = '%Y %m %d %H %M'
    col_offset = 0
    date_length = 16
  elif  obs[0].find('YY  MM DD hh mm') >= 0:
    obs_frmt = '%Y %m %d %H %M'
    col_offset = 0
    date_length = 16
  elif  obs[0].find('YYYY MM DD hh WD') >= 0:
    obs_frmt = '%Y %m %d %H'
    col_offset = -1
    date_length = 13
  else:
    print( 'problem with observed data format')
    raise SystemExit(0)
   
  if output_datetime is not None:
    # Get data at specified datetime values
    nlines = len(obs)
    lines_searched = 0
    for time in output_datetime:
      found = False
      obs_data['datetime'].append(time)
      for j in range(lines_searched,nlines):
        line = obs[j]
        if line.find('#') >= 0 or len(line.strip()) == 0 or not line[0].isdigit():
          continue
        date = line[0:date_length]
        date_time = datetime.datetime.strptime(date,obs_frmt)
        if date_time == time:
          for var in variables:
            col = variables[var]['obs_col'] + col_offset
            obs_data[var].append(line.split()[col])
          lines_searched = j+1
          found = True
          break
        if date_time > time:
          lines_searched = j
          break
        if j == nlines-1 and date_time < time:
          lines_searched = j
          break
      if found == False:
        for var in variables:
          obs_data[var].append(variables[var]['fill_val'])
      #print(obs_data['datetime'][-1],obs_data['hs'][-1])
  else:
    # Get data from observation file between min and max output times
    for line in obs[1:]:
      if line.find('#') >= 0 or len(line.strip()) == 0 or not line[0].isdigit():
        continue
      date = line[0:date_length]
      date_time = datetime.datetime.strptime(date,obs_frmt)
      if date_time >= datetime.datetime.strptime(min_date,frmt) and \
         date_time <= datetime.datetime.strptime(max_date,frmt):
        obs_data['datetime'].append(date_time)
        for var in variables:
          col = variables[var]['obs_col'] + col_offset
          obs_data[var].append(line.split()[col])

  # Convert observation data and replace fill values with nan
  for var in variables:
    obs_data[var] = np.asarray(obs_data[var])
    obs_data[var] = obs_data[var].astype(float)
    fill_val = variables[var]['fill_val']
    obs_data[var][obs_data[var] >= fill_val] = np.nan

  obs_data['datetime'] = np.asarray(obs_data['datetime'],dtype='O')

  return obs_data

################################################################################################
################################################################################################

def output_time_to_date(output_time,ref_date,start_date=None):

  output_date = []
  output_datetime = []

  date_min = '1000 01 01 00 00'
  frmt = '%Y %m %d %H %M'
  adjust_start = False
  if ref_date + datetime.timedelta(days=output_time[0]) < datetime.datetime.strptime(date_min,frmt):
    adjust_start = True

  # Handle reference time
  if start_date and adjust_start:                                             # start_time is used to create a common starting
    start_date = datetime.datetime.strptime(start_date,'%Y-%m-%d %H:%M:%S')     # point across runs that start from different initial
    offset = -datetime.timedelta(days=output_time[0])                                  # dates i.e. E3SM 0001-01-01, WW3 2005-06-01
  else:
    start_date = ref_date
    offset = datetime.timedelta(days=0)

  # Convert output times to date format
  for i,t in enumerate(output_time):
    try:
      date = start_date + datetime.timedelta(days=t) + offset
    except:
      print( "likely datetime overflow error")
      print( "  try specifying start_date in the config file")
      print( "  if simulation started at something like 0001-01-01")
      raise SystemExit(0)
     
    output_date.append(date.strftime('%Y %m %d %H %M'))
    output_datetime.append(date)
  output_datetime = np.asarray(output_datetime,dtype='O')

  return output_datetime, output_date

################################################################################################
################################################################################################

if __name__ == '__main__':

  pwd = os.getcwd()
  
  # Read config file  
  f = open(pwd+'/plot_points.config')
  cfg = yaml.load(f,Loader=yaml.Loader)
  pprint.pprint(cfg)

  if 'start_date' not in cfg:
    cfg['start_date'] = None
  
  # Create list of output files for each run
  runs = {}
  field = {}
  for run in cfg['model_direcs']:
    direc = cfg['model_direcs'][run]
    if direc[-1] != '/':
      direc = direc + '/'
    if not os.path.exists(direc):
      print('Directory for {} run data not found'.format(run))
      raise SystemExit(0)
    wav_files = sorted(glob.glob(direc+'ww3*_tab.nc'))
    wnd_files = sorted(glob.glob(direc+'cfsr*_tab.nc'))
    if len(wav_files) > 0:
      runs[run] = [wav_files,wnd_files]
      field[run] = False
    else:
      wav_files = sorted(glob.glob(direc+'*ww3*.nc'))
      runs[run] = [wav_files]
      field[run] = True
  
  
  #-----------------------------------
  # Read point data from NetCDF files
  #-----------------------------------
  data = {}
  stations = {}
  if len(runs) > 0:
    for k,run in enumerate(runs):
      data_files = runs[run]
      if field[run]:
        data[run],stations[run] = interpolate_stations_from_fields(data_files,variables,cfg["station_file"],cfg['start_date']) 
      else:
        data[run],stations[run] = read_point_files(data_files,variables,cfg['start_date'])  
    
    # Get list of all stations
    station_list = [] 
    for run in stations:
      for sta in stations[run]['name']:
        station_list.append(sta)
    station_list = list(set(station_list))

    # Find overall date range of data
    date_min = '3000 01 01 00 00'
    date_max = '1000 01 01 00 00'
    frmt = '%Y %m %d %H %M'
    for run in data:
      if data[run]['datetime'][0]  < datetime.datetime.strptime(date_min,frmt):
        date_min = data[run]['date'][0]
      if data[run]['datetime'][-1] > datetime.datetime.strptime(date_max,frmt):
        date_max = data[run]['date'][-1]
    print( date_min,date_max)
  else:
    stations["stations_only"] = read_station_file(cfg["station_file"])
    station_list = []
    for sta in stations["stations_only"]["name"]:
      station_list.append(sta) 
   
    date_min = cfg["date_min"] 
    date_max = cfg["date_max"] 
    
  #-------------------------------------------------
  # Read observation data and plot for each station
  #-------------------------------------------------
  for i,sta in enumerate(station_list):
    print( sta)
  
    # Check if observation file exists
    obs_file = ""
    obs_file_check = cfg['obs_direc']+sta+'_'+cfg['year']+'.txt'
    if os.path.isfile(obs_file_check):
      obs_file = obs_file_check
    
    obs_file_check = cfg['obs_direc']+sta+'_'+cfg['year']+'_stdmet.txt'
    if os.path.isfile(obs_file_check):
      obs_file = obs_file_check

    obs_file_check = cfg['obs_direc']+sta+'.txt'
    if os.path.isfile(obs_file_check):
      obs_file = obs_file_check


    # Get data from observation file at output times
    if not obs_file:      
      print( '  no observation file found')
    else:
      print( '  '+obs_file)
      obs_data = read_station_data(obs_file,date_min,date_max,variables)
        
      # Create figure 
      fig = plt.figure(figsize=[6,2+2*len(variables)])
      gs = gridspec.GridSpec(nrows=len(variables)+1,ncols=2,figure=fig)
     
      # Find station location 
      for run in stations:
        if sta in stations[run]['name']:
          ind = stations[run]['name'].index(sta)
          lon = stations[run]['lon'][ind]
          lat = stations[run]['lat'][ind]
          break
     
      # Plot global station location
      ax = fig.add_subplot(gs[0,0],projection=ccrs.PlateCarree())
      ax.set_extent([-180.0, 180.0, -90.0, 90.0],crs=ccrs.PlateCarree())
      ax.add_feature(cfeature.LAND, zorder=100)
      ax.add_feature(cfeature.LAKES, alpha=0.5, zorder=101)
      ax.coastlines(zorder=101)
      ax.plot(lon,lat,'ro',zorder=102)
      
      # Plot local station location
      ax = fig.add_subplot(gs[0,1],projection=ccrs.PlateCarree())
      ax.set_extent([lon-10.0, lon+10.0, lat-7.0, lat+7.0],crs=ccrs.PlateCarree())
      ax.add_feature(cfeature.LAND, zorder=100)
      ax.add_feature(cfeature.LAKES, alpha=0.5, zorder=101)
      ax.coastlines(zorder=101)
      ax.plot(lon,lat,'ro',zorder=102)
  
      # Plot modeled and observed data timeseries
      for k,var in enumerate(variables):
        print( '  '+var)
        lines = []
        labels = []
        data_plotted = False
        ax = fig.add_subplot(gs[k+1,:])
        if np.size(obs_data['datetime']) > 0:
          l1, = ax.plot(obs_data['datetime'],obs_data[var])
          lines.append(l1)
          labels.append('Observed')
          data_plotted = True
        for run in data:
          if sta in stations[run]['name']:
            ind = stations[run]['name'].index(sta)
            if variables[var]['recip'] == True:
              data[run][var][:,ind] = 1.0/data[run][var][:,ind]
            l2, = ax.plot(data[run]['datetime'],data[run][var][:,ind])
            lines.append(l2)
            labels.append(run)
            data_plotted = True
        ax.set_title(variables[var]['title'])
        if data_plotted:
          ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
        ax.xaxis.set_major_locator(mdates.WeekdayLocator(interval=1))
        ax.set_xlabel('time')
        ax.set_ylabel(variables[var]['label']+' ('+variables[var]['units']+')')
        #print ax.get_xlim()
        if variables[var]['units'] == 'deg':
          ax.set_ylim([0.0,360.0])          
        if ax.get_xlim() == (-0.001, 0.001): # Detect when there is no availiable data
          fig.delaxes(ax)                    # And delete axis to prevent an error
      lgd = plt.legend(lines,labels,loc=9,bbox_to_anchor=(0.5,-0.5),ncol=2,fancybox=False,edgecolor='k')
      st = plt.suptitle('Station '+sta,y=1.025,fontsize=16)
      fig.canvas.draw()
      fig.tight_layout()
      fig.savefig(sta+'.png',bbox_inches='tight',bbox_extra_artists=(lgd,st,))
      plt.close()
  
  if not os.path.exists(cfg["plot_direc"]):
    subprocess.call(['mkdir','-p',cfg["plot_direc"]])
  subprocess.call('mv *.png '+cfg["plot_direc"],shell=True)
