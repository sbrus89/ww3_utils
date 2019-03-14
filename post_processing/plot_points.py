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
from mpl_toolkits.basemap import Basemap
from scipy import interpolate
plt.switch_backend('agg')
np.set_printoptions(threshold=np.nan)
#--------------------------
# Define variables to plot
#--------------------------
variables = {'hs'  :{'obs_col' : 8,
                     'fill_val': 99.00,
                     'recip'   : False,
                     'label'   : 'Significant wave height',
                     'units'   : 'm',
                     'aka'     : ['HS']},
             'th1p':{'obs_col' : 11, 
                     'fill_val': 999,
                     'recip'   : False,
                     'label'   : 'Dominant wave direction',
                     'units'   : 'deg',
                     'aka'     : ['dp']},
             'fp'  :{'obs_col' : 9,
                     'fill_val': 99.0,
                     'recip'   : True,
                     'label'   : 'Dominant wave period',
                     'units'   : 's',
                     'aka'     : ['FP']},
             'wnd' :{'obs_col' : 6,
                     'fill_val': 99.0,
                     'recip'   : False,
                     'label'   : 'Wind speed',
                     'units'   : 'm/s',
                     'aka'     : ['uwnd','vwnd']},
             'wnddir':{'obs_col' : 5,
                       'fill_val': 999,
                       'recip'   : False,
                       'label'   : 'Wind direction',
                       'units'   : 'deg',
                       'aka'     : ['uwnd','vwnd']}}

################################################################################################
################################################################################################

def read_point_files(data_files,variables):

  for j,files in enumerate(data_files):
    for i,name in enumerate(files):
      print name.split('/')[-1]
    
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
          print var
          data[var][i][:] = nc_file.variables[var][:]

  data['datetime'], data['date'] = output_time_to_date(data['time'],data['ref_date'])

  return data, stations

################################################################################################
################################################################################################

def interpolate_stations_from_fields(data_files,variables,station_file):

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
  sta_pts = np.column_stack((stations['lon'],stations['lat']))

  t = 0.0
  for j,files in enumerate(data_files):
    for i,name in enumerate(files):
      print name.split('/')[-1]

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
      Lon,Lat = np.meshgrid(lon,lat)

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
          print var      

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

  data['datetime'], data['date'] = output_time_to_date(data['time'],data['ref_date'])

  return data, stations
################################################################################################
################################################################################################

def read_field(nc_file,var):
  
  if 'time' in nc_file.variables:
    field = nc_file.variables[var][0,:,:]
  else:
    field = nc_file.variables[var][:,:]

  return field

################################################################################################
################################################################################################

def read_station_data(obs_file,min_date,max_date,variables):

  frmt = '%Y %m %d %H %M'

  # Initialize variable for observation data
  obs_data = {}
  for var in variables:
    obs_data[var] = []  
  obs_data['datetime'] = []

  # Get data from observation file between min and max output times
  f = open(obs_file)
  obs = f.read().splitlines()
  for line in obs[1:]:
    date = line[0:16]
    date_time = datetime.datetime.strptime(date,frmt)
    if date_time >= datetime.datetime.strptime(min_date,frmt) and \
       date_time <= datetime.datetime.strptime(max_date,frmt):
      obs_data['datetime'].append(date_time)
      for var in variables:
        col = variables[var]['obs_col']
        obs_data[var].append(line.split()[col])

  # Convert observation data and replace fill values with nan
  for var in variables:
    obs_data[var] = np.asarray(obs_data[var])
    obs_data[var] = obs_data[var].astype(np.float)
    fill_val = variables[var]['fill_val']
    obs_data[var][obs_data[var] >= fill_val] = np.nan

  obs_data['datetime'] = np.asarray(obs_data['datetime'],dtype='O')

  return obs_data

################################################################################################
################################################################################################

def output_time_to_date(output_time,ref_date):

  output_date = []
  output_datetime = []

  # Handle reference time
  if 'start_date' in cfg:                                                              # start_time is used to create a common starting
    start_date = datetime.datetime.strptime(cfg['start_date'],'%Y-%m-%d %H:%M:%S')     # point across runs that start from different initial
    offset = -datetime.timedelta(days=output_time[0])                                  # dates i.e. E3SM 0001-01-01, WW3 2005-06-01
  else:
    start_date = ref_date
    offset = datetime.timedelta(days=0)

  # Convert output times to date format
  for i,t in enumerate(output_time):
    try:
      date = start_date + datetime.timedelta(days=t) + offset
    except:
      print "likely datetime overflow error"
      print "  try specifying start_date in the config file"
      print "  if simulation started at something like 0001-01-01"
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
  cfg = yaml.load(f)
  pprint.pprint(cfg)
  
  # Create list of output files for each run
  runs = {}
  field = {}
  for run in cfg['model_direcs']:
    direc = cfg['model_direcs'][run]
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
  for k,run in enumerate(runs):
    data_files = runs[run]
    if field[run]:
      data[run],stations[run] = interpolate_stations_from_fields(data_files,variables,cfg["station_file"]) 
    else:
      data[run],stations[run] = read_point_files(data_files,variables)  
  
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
  print date_min,date_max
    
    
  #-------------------------------------------------
  # Read observation data and plot for each station
  #-------------------------------------------------
  for i,sta in enumerate(station_list):
    print sta
  
  
    # Get data from observation file at output times
    obs_file = cfg['obs_direc']+sta+'_'+cfg['year']+'.txt'
    if os.path.isfile(obs_file):
      obs_data = read_station_data(obs_file,date_min,date_max,variables)
        
      # Create figure 
      fig = plt.figure(figsize=[6,12])
      gs = gridspec.GridSpec(nrows=len(variables)+1,ncols=2,figure=fig)
     
      # Find station location 
      for run in stations:
        if sta in stations[run]['name']:
          ind = stations[run]['name'].index(sta)
          lon = stations[run]['lon'][ind]
          lat = stations[run]['lat'][ind]
          break
     
      # Plot global station location
      ax = fig.add_subplot(gs[0,0])
      m = Basemap(projection='cyl',llcrnrlat= -90,urcrnrlat=90,\
                                   llcrnrlon=-180,urcrnrlon=180,resolution='c')
      m.fillcontinents(color='tan',lake_color='lightblue')
      m.drawcoastlines()
      ax.plot(lon,lat,'ro')
      
      # Plot local station location
      ax = fig.add_subplot(gs[0,1])
      m = Basemap(projection='cyl',llcrnrlat=lat-7.0 ,urcrnrlat=lat+7.0,\
                                   llcrnrlon=lon-10.0,urcrnrlon=lon+10.0,resolution='l')
      m.fillcontinents(color='tan',lake_color='lightblue')
      m.drawcoastlines()
      ax.plot(lon,lat,'ro')
  
      # Plot modeled and observed data timeseries
      for k,var in enumerate(variables):
        print '  '+var
        lines = []
        labels = []
        ax = fig.add_subplot(gs[k+1,:])
        l1, = ax.plot(obs_data['datetime'],obs_data[var])
        lines.append(l1)
        labels.append('Observed')
        for run in data:
          if sta in stations[run]['name']:
            ind = stations[run]['name'].index(sta)
            if variables[var]['recip'] == True:
              data[run][var][:,ind] = 1.0/data[run][var][:,ind]
            l2, = ax.plot(data[run]['datetime'],data[run][var][:,ind])
            lines.append(l2)
            labels.append(run)
        ax.set_title(variables[var]['label'])
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
        #ax.xaxis.set_major_locator(plt.MaxNLocator(6))
        ax.set_xlabel('time')
        ax.set_ylabel(var)
        if variables[var]['units'] == 'deg':
          ax.set_ylim([0.0,360.0])          
        if ax.get_xlim() == (-0.001, 0.001): # Detect when there is no availiable data
          fig.delaxes(ax)                    # And delete axis to prevent an error
      lgd = plt.legend(lines,labels,loc=9,bbox_to_anchor=(0.5,-0.5),ncol=2,fancybox=False,edgecolor='k')
      st = plt.suptitle('Station '+sta,y=1.025,fontsize=16)
      fig.tight_layout()
      fig.savefig(sta+'.png',bbox_inches='tight',bbox_extra_artists=(lgd,st,))
      plt.close()
  
  if not os.path.exists(cfg["plot_direc"]):
    subprocess.call(['mkdir','-p',cfg["plot_direc"]])
  subprocess.call('mv *.png '+cfg["plot_direc"],shell=True)
