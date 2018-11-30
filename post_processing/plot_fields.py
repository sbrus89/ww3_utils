import netCDF4 
import matplotlib.pyplot as plt

from mpl_toolkits.basemap import Basemap
import numpy as np
import glob
import datetime
import os
import yaml
import pprint
import subprocess
from mpl_toolkits.basemap import Basemap
plt.switch_backend('agg')

def read_field_ww3(f):

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

  var = nc_file.variables['hs'][0,:,:]

  return lon,lat,var,output_date


def read_field_e3sm(f):

  nc_file = netCDF4.Dataset(f,'r')

  output_date = f.split('.')[-2] 

  lon  = np.linspace(cfg["lon_range"][0],cfg["lon_range"][1],nc_file.dimensions['NX'].size)
  lat  = np.linspace(cfg["lat_range"][0],cfg["lat_range"][1],nc_file.dimensions['NY'].size)

  var = nc_file.variables['hs'][:,:]

  return lon,lat,var,output_date

###############################################################################################
###############################################################################################

pwd = os.getcwd()


# Read in configuration file
f = open(pwd+'/plot_fields.config')
cfg = yaml.load(f)
pprint.pprint(cfg)

# Setup to read either e3sm or ww3 output
if 'model' in cfg and cfg['model'] == 'e3sm':
  read_field = read_field_e3sm
  file_wildcard = '*.ww3.*.nc'
else:
  read_field = read_field_ww3
  file_wildcard = '*.nc'

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
elif nruns == 2:
  cmap = 'bwr'
  diff = ' difference'

# Get field output files
output_direc = cfg['model_direc'][0][1]
files = sorted(glob.glob(output_direc+file_wildcard))

count = 0
for f in files:
  filename = f.split('/')[-1]
  print filename

  # Read in and compute fields to plot
  lon,lat,var,output_date = read_field(f)
  if nruns > 1:
    lon2,lat2,var2,output_date2 = read_field(cfg['model_direc'][1][1]+filename)
    if var.shape == var2.shape:
      var = var-var2
    elif (var.shape[0] == var2.shape[0]/2) and (var.shape[1] == var2.shape[1]/2):
      var = var-var2[::2,::2]
    elif (var.shape[0] == var2.shape[0]*2) and (var.shape[1] == var2.shape[1]*2):
      var = var[::2,::2]-var2
      lon = lon2
      lat = lat2

  # Accumulate for average
  count = count + 1
  if count == 1:
    var_avg = np.copy(var)
  else:
    var_avg = var_avg + var
  
  # Plot timesnap
  if cfg["timesnaps"]:    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    m = Basemap(projection='cyl',llcrnrlat= -90,urcrnrlat=90,\
                                 llcrnrlon=0,urcrnrlon=360,resolution='c')
    m.fillcontinents(color='tan',lake_color='lightblue')
    m.drawcoastlines()
    if "tsnap_range" in cfg:
      levels = np.linspace(cfg["tsnap_range"][0],cfg["tsnap_range"][1],100)
      cf = ax.contourf(lon,lat,var,levels=levels,cmap=cmap)
    else: 
      cf = ax.contourf(lon,lat,var,100,cmap=cmap)
    ax.set_title('Wave height'+diff+' on '+output_date)
    cbar = plt.colorbar(cf,ax=ax,orientation='horizontal')
    cbar.set_label('wave height'+diff+' (m)')
    plt.savefig('-'.join(output_date.split())+'.png')
    plt.close() 
  
# Compute average
var_avg = var_avg/float(count)

# Plot average
fig = plt.figure()
m = Basemap(projection='cyl',llcrnrlat= -90,urcrnrlat=90,\
                             llcrnrlon=0,urcrnrlon=360,resolution='c')
m.fillcontinents(color='tan',lake_color='lightblue')
m.drawcoastlines()
if "avg_range" in cfg:
  levels = np.linspace(cfg["avg_range"][0],cfg["avg_range"][1],100)
  cf = plt.contourf(lon,lat,var_avg,levels=levels,cmap=cmap)
else:
  cf = plt.contourf(lon,lat,var_avg,100,cmap=cmap)
plt.title('Average wave height'+diff)
cbar = plt.colorbar(cf,orientation='horizontal')
cbar.set_label('wave height'+diff+' (m)')
plt.savefig('average.png')
plt.close() 

# Move plots to their own directory
if not os.path.exists(cfg["plot_direc"]):
  subprocess.call(['mkdir','-p',cfg["plot_direc"]])
subprocess.call('mv *.png '+cfg["plot_direc"],shell=True)    
