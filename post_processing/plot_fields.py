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

pwd = os.getcwd()

# Read in configuration file
f = open(pwd+'/plot_fields.config')
cfg = yaml.load(f)
pprint.pprint(cfg)

# Find nc file names
output_direc = cfg['model_direc']['direc']
files = sorted(glob.glob(output_direc+'*.nc'))

for f in files:
  print f.split('/')[-1]

  nc_file = netCDF4.Dataset(f,'r')

  ref_date = nc_file.variables['time'].getncattr('units').replace('days since ','')
  ref_date = datetime.datetime.strptime(ref_date,'%Y-%m-%d %H:%M:%S')
  output_time = nc_file.variables['time'][:]

  output_date = []
  for t in output_time:
    date = ref_date + datetime.timedelta(days=t)
    output_date.append(date.strftime('%Y %m %d %H %M'))

  lon  = nc_file.variables['longitude'][:] 
  lat  = nc_file.variables['latitude'][:]

  for i,t in enumerate(output_date): 
    var  = nc_file.variables['hs'][i,:,:]
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    levels = np.linspace(0.0,10.0,100)
    cf = ax.contourf(lon,lat,var,levels=levels)
    ax.axis('image')
    ax.set_title(t)
    plt.colorbar(cf,ax=ax,orientation='horizontal')
    plt.savefig('-'.join(t.split())+'.png')
    plt.close() 
  
if not os.path.exists(cfg["plot_direc"]):
  subprocess.call(['mkdir','-p',cfg["plot_direc"]])
  subprocess.call('mv *.png '+cfg["plot_direc"],shell=True)    
