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
plt.switch_backend('agg')

#--------------------------
# Define variables to plot
#--------------------------
variables = {'hs'  :{'obs_col' : 8,
                     'fill_val': 99.00,
                     'recip'   : False,
                     'label'   : 'Significant wave height'},
             'th1p':{'obs_col' : 11, 
                     'fill_val': 999,
                     'recip'   : False,
                     'label'   : 'Dominant wave direction'},
             'fp'  :{'obs_col' : 9,
                     'fill_val': 99.0,
                     'recip'   : True,
                     'label'   : 'Dominant wave period'},
             'wnd' :{'obs_col' : 6,
                     'fill_val': 99.0,
                     'recip'   : False,
                     'label'   : 'Wind speed'},
             'wnddir':{'obs_col' : 5,
                       'fill_val': 999,
                       'recip'   : False,
                       'label'   : 'Wind direction'}}

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
          output_time = np.empty(nfiles)
          ref_date = nc_file.variables['time'].getncattr('units').replace('days since ','')
          ref_date = datetime.datetime.strptime(ref_date,'%Y-%m-%d %H:%M:%S')
      
      # Get time and output variables
      if j == 0:
        output_time[i] = nc_file.variables['time'][:]
      for var in variables:
        if var in nc_file.variables:
          data[var][i][:] = nc_file.variables[var][:]

  return data, stations, output_time, ref_date

################################################################################################
################################################################################################

def read_station_data(obs_file,output_date,variables):

  # Initialize variable for observation data
  obs_data = {}
  for var in variables:
    obs_data[var] = []  

  # Get data from observation file at output times
  f = open(obs_file)
  obs = f.read().splitlines()
  nlines = len(obs)
  lines_searched = 0
  for t in output_date:
    found = False
    for j in range(lines_searched,nlines):
      if obs[j].find(t) >= 0:
        for var in variables:
          col = variables[var]['obs_col']
          obs_data[var].append(obs[j].split()[col])
        lines_searched = j+1
        found = True
        break
    if found == False:
      for var in variables:
        obs_data[var].append('999.0')

  # Convert observation data and replace fill values with nan
  for var in variables:
    obs_data[var] = np.asarray(obs_data[var])
    obs_data[var] = obs_data[var].astype(np.float)
    fill_val = variables[var]['fill_val']
    obs_data[var][obs_data[var] >= fill_val] = np.nan

  return obs_data

################################################################################################
################################################################################################

if __name__ == '__main__':

  pwd = os.getcwd()
  
  
  f = open(pwd+'/plot_points.config')
  cfg = yaml.load(f)
  pprint.pprint(cfg)
  
  
  runs = {}
  for run in cfg['model_direcs']:
    direc = cfg['model_direcs'][run]
    wav_files = sorted(glob.glob(direc+'ww3*_tab.nc'))
    wnd_files = sorted(glob.glob(direc+'cfsr*_tab.nc'))
    runs[run] = [wav_files,wnd_files]
  
  
  
  #-----------------------------------
  # Read point data from NetCDF files
  #-----------------------------------
  data = {}
  stations = {}
  for k,run in enumerate(runs):
    data_files = runs[run]
    data[run],stations[run],output_time,ref_date = read_point_files(data_files,variables)  
  
  # Convert output times to date format
  output_date = []
  output_datetime = []
  for t in output_time:
    date = ref_date + datetime.timedelta(days=t)
    output_date.append(date.strftime('%Y %m %d %H %M'))
    output_datetime.append(date)
  output_datetime = np.asarray(output_datetime,dtype='O')
  
  #-------------------------------------------------
  # Read observation data and plot for each station
  #-------------------------------------------------
  station_list = [] 
  for run in stations:
    for sta in stations[run]['name']:
      station_list.append(sta)
  station_list = list(set(station_list))
  
  for i,sta in enumerate(station_list):
    print sta
  
  
    # Get data from observation file at output times
    obs_file = cfg['obs_direc']+sta+'_'+cfg['year']+'.txt'
    if os.path.isfile(obs_file):
      obs_data = read_station_data(obs_file,output_date,variables)
        
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
  
      # Determine if model data is not available, i.e. values are all the same  (this causes issues with the plot axis)    
      count = 0
      for run in data:
        flag = False
        for var in variables:
          if sta in stations[run]['name']:
            ind = stations[run]['name'].index(sta)
            if len(np.unique(data[run][var][:,ind])) == 1:
              flag = True
        if flag == True:
          count = count + 1

      plot_flag = True
      if count == len(data):
        plot_flag = False
  
      if plot_flag == False:
  
        print "  model data not availiable"
        st = plt.suptitle('Station '+sta,y=1.025,fontsize=16)
        fig.tight_layout()
        fig.savefig(sta+'.png',bbox_inches='tight')
        plt.close()
        continue
  
      else:
  
        # Plot modeled and observed data timeseries
        for k,var in enumerate(variables):
          print '  '+var
          lines = []
          labels = []
          ax = fig.add_subplot(gs[k+1,:])
          l1, = ax.plot(output_datetime,obs_data[var])
          lines.append(l1)
          labels.append('Observed')
          for run in data:
            if sta in stations[run]['name']:
              ind = stations[run]['name'].index(sta)
              if variables[var]['recip'] == True:
                data[run][var][:,ind] = 1.0/data[run][var][:,ind]
              l2, = ax.plot(output_datetime,data[run][var][:,ind])
              lines.append(l2)
              labels.append(run)
          ax.set_title(variables[var]['label'])
          ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
          #ax.xaxis.set_major_locator(plt.MaxNLocator(6))
          ax.set_xlabel('time')
          ax.set_ylabel(var)
        lgd = plt.legend(lines,labels,loc=9,bbox_to_anchor=(0.5,-0.5),ncol=2,fancybox=False,edgecolor='k')
        st = plt.suptitle('Station '+sta,y=1.025,fontsize=16)
        fig.tight_layout()
        fig.savefig(sta+'.png',bbox_inches='tight',bbox_extra_artists=(lgd,st,))
        plt.close()
  
  if not os.path.exists(cfg["plot_direc"]):
    subprocess.call(['mkdir','-p',cfg["plot_direc"]])
  subprocess.call('mv *.png '+cfg["plot_direc"],shell=True)
