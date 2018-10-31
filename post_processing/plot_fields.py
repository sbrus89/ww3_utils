import netCDF4 
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import glob
import datetime
import os
import yaml
import pprint
import subprocess
from mpl_toolkits.basemap import Basemap
plt.switch_backend('agg')

def read_field(f):

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

###############################################################################################
###############################################################################################

pwd = os.getcwd()

# Read in configuration file
f = open(pwd+'/plot_fields.config')
cfg = yaml.load(f)
pprint.pprint(cfg)

# Find nc file names
runs = [x for x in cfg['model_direc']]
nruns = len(runs)
if nruns > 2:
  print "Only 2 runs can be specified for difference plots"
  raise SystemExit(0)
output_direc = cfg['model_direc'][runs[0]]
files = sorted(glob.glob(output_direc+'*.nc'))

count = 0
for f in files:
  filename = f.split('/')[-1]
  print filename

  lon,lat,var,output_date = read_field(f)
  if nruns > 1:
    lon2,lat2,var2,output_date2 = read_field(cfg['model_direc'][runs[1]]+filename)
    if var.shape == var2.shape:
      var = np.absolute(var-var2)
    elif (var.shape[0] == var2.shape[0]/2) and (var.shape[1] == var2.shape[1]/2):
      var = np.absolute(var-var2[::2,::2])
    elif (var.shape[0] == var2.shape[0]*2) and (var.shape[1] == var2.shape[1]*2):
      var = np.absolute(var[::2,::2]-var2)
      lon = lon2
      lat = lat2

  count = count + 1
  if count == 1:
    var_avg = np.copy(var)
  else:
    var_avg = var_avg + var
      
  fig = plt.figure()
  ax = fig.add_subplot(1,1,1)
  cf = ax.contourf(lon,lat,var,100)
  ax.axis('image')
  ax.set_title(output_date)
  plt.colorbar(cf,ax=ax,orientation='horizontal')
  plt.savefig('-'.join(output_date.split())+'.png')
  plt.close() 
  
var_avg = var_avg/float(count)
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
cf = ax.contourf(lon,lat,var_avg,100)
ax.axis('image')
ax.set_title('Average')
plt.colorbar(cf,ax=ax,orientation='horizontal')
plt.savefig('average.png')
plt.close() 

if not os.path.exists(cfg["plot_direc"]):
  subprocess.call(['mkdir','-p',cfg["plot_direc"]])
subprocess.call('mv *.png '+cfg["plot_direc"],shell=True)    
