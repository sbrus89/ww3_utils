import netCDF4 
import matplotlib.pyplot as plt
import numpy as np
import glob
import pprint
import datetime
import os
plt.switch_backend('agg')

files = sorted(glob.glob('*.nc'))
nfiles = len(files)
year = '2005'


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
                     'label'   : 'Dominant wave period'}}

for i,name in enumerate(files):
  print name

  nc_file = netCDF4. Dataset(name,'r')

  # Initializations for first iteration
  if i == 0:
    stations = netCDF4.chartostring(nc_file.variables['station_name'][:,:])
    nstations = len(stations)
    
    data = {}
    for var in variables:
      data[var] = np.empty((nfiles,nstations))

    output_time = np.empty(nfiles)
    ref_date = nc_file.variables['time'].getncattr('units').replace('days since ','')
    ref_date = datetime.datetime.strptime(ref_date,'%Y-%m-%d %H:%M:%S')
  
  # Get time and output variables
  output_time[i] = nc_file.variables['time'][:]
  for var in variables:
    data[var][i][:] = nc_file.variables[var][:]

# Convert output times to date format
output_datetime = []
for t in output_time:
  date = ref_date + datetime.timedelta(days=t)
  output_datetime.append(date.strftime('%Y %m %d %H %M'))


for i,sta in enumerate(stations):
  print sta

  # Initialize variable for observation data
  obs_data = {}
  for var in variables:
    obs_data[var] = []  

  # Get data from observation file at output times
  obs_file = '../ww3_utils/stations/'+sta+'_'+year+'.txt'
  if os.path.isfile(obs_file):
    f = open(obs_file)
    obs = f.read().splitlines()
    nlines = len(obs)
    lines_searched = 0
    for t in output_datetime:
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

    # Convert observation data
    for var in variables:
      obs_data[var] = np.asarray(obs_data[var])
      obs_data[var] = obs_data[var].astype(np.float)
      plot_obs_data = False
      fill_val = variables[var]['fill_val']
      if np.amin(obs_data[var]) < fill_val: 
        plot_obs_data = True
      obs_data[var][obs_data[var] >= fill_val] = np.nan
      if variables[var]['recip'] == True:
        obs_data[var] = 1.0/obs_data[var]

    plt.figure()
    for k,var in enumerate(variables):
      print '  '+var
      plt.subplot(len(variables),1,k+1)
      plt.plot(output_time,data[var][:,i])
      plt.plot(output_time,obs_data[var])

      #plt.title('Station '+sta)
      #plt.xlabel('time')
      #plt.ylabel(var)
    plt.savefig(sta+'_'+var+'.png',bbox_inches='tight')
    plt.close()


