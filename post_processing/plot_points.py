import netCDF4 
import matplotlib.pyplot as plt
import numpy as np
import glob
import pprint
import datetime
import os
from mpl_toolkits.basemap import Basemap
plt.switch_backend('agg')

pwd = os.getcwd()

model_direc = pwd+'/model_data/'
obs_direc = pwd+'/obs_data/'

wav_files = sorted(glob.glob(model_direc+'ww3*.nc'))
wnd_files = sorted(glob.glob(model_direc+'cfsr*.nc'))
nfiles = len(wav_files)
year = '2005'

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


#-----------------------------------
# Read point data from NetCDF files
#-----------------------------------
for j,files in enumerate([wav_files,wnd_files]):
  for i,name in enumerate(files):
    print name.split('/')[-1]
  
    nc_file = netCDF4.Dataset(name,'r')
  
    # Initializations for first iteration
    if j == 0 and i == 0:
        # Station names
        stations = netCDF4.chartostring(nc_file.variables['station_name'][:,:])
        nstations = len(stations)
      
        # Station locations
        station_loc = {}
        station_loc['lon'] = np.squeeze(nc_file.variables['longitude'][:])
        station_loc['lat'] = np.squeeze(nc_file.variables['latitude'][:])
  
        # Model output data
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


# Convert output times to date format
output_date = []
output_datetime = []
for t in output_time:
  date = ref_date + datetime.timedelta(days=t)
  output_date.append(date.strftime('%Y %m %d %H %M'))
  output_datetime.append(date)

#-------------------------------------------------
# Read observation data and plot for each station
#-------------------------------------------------
for i,sta in enumerate(stations):
  print sta

  # Initialize variable for observation data
  obs_data = {}
  for var in variables:
    obs_data[var] = []  

  # Get data from observation file at output times
  obs_file = obs_direc+sta+'_'+year+'.txt'
  if os.path.isfile(obs_file):
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
      plot_obs_data = False
      fill_val = variables[var]['fill_val']
      if np.amin(obs_data[var]) < fill_val: 
        plot_obs_data = True
      obs_data[var][obs_data[var] >= fill_val] = np.nan

    # Plot station location
    fig = plt.figure(figsize=[6,12])
    ax = fig.add_subplot(len(variables)+1,1,1)
    m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
                llcrnrlon=-180,urcrnrlon=180,resolution='c')
    m.fillcontinents(color='tan',lake_color='lightblue')
    m.drawcoastlines()
    ax.plot(station_loc['lon'][i],station_loc['lat'][i],'ro')
    ax.set_title('Station '+sta)

    # Plot modeled and observed data timeseries
    for k,var in enumerate(variables):
      print '  '+var
      ax = fig.add_subplot(len(variables)+1,1,k+2)
      if variables[var]['recip'] == True:
        data[var][:,i] = 1.0/data[var][:,i]
      l1, = ax.plot(output_datetime,data[var][:,i])
      l2, = ax.plot(output_datetime,obs_data[var])
      ax.set_title(variables[var]['label'])
      ax.set_xlabel('time')
      ax.set_ylabel(var)
    lgd = plt.legend((l1,l2),('Modeled','Observed'),loc=9,bbox_to_anchor=(0.5,-0.5),ncol=2,fancybox=False,edgecolor='k')
    fig.tight_layout()
    fig.savefig(sta+'_'+var+'.png',bbox_inches='tight',bbox_extra_artists=(lgd,))
    plt.close()


