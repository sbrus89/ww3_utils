import plot_points
from scipy import stats
from scipy import interpolate
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.cm as cm
import matplotlib.colors as colors
import glob
import datetime
import numpy as np
import os
import netCDF4 as nc4
from collections import OrderedDict

direc = OrderedDict()
direc['1/2 degree'] = '/users/sbrus/scratch4/WW3_testing/glo_30m/post_processing/june-august/model_data/points/'
direc['1 degree']   = '/users/sbrus/scratch4/WW3_testing/glo_1d/post_processing/june-august/model_data/points/'
direc['2 degree']   = '/users/sbrus/scratch4/WW3_testing/glo_2d/post_processing/june-august/model_data/points/'

obs_direc = './obs_data/'
variables = plot_points.variables
variable = 'hs'
year = '2005'
bathy_file = '/users/sbrus/climate/bathy_data/SRTM15_plus/earth_relief_15s.nc'

# Read in modeled data
data = {}
stations = {}
for run in direc:
  wav_files = [sorted(glob.glob(direc[run]+'ww3*_tab.nc'))]
  data[run], stations[run], output_time, ref_date = plot_points.read_point_files(wav_files,variables)

# Convert output times to date format
output_date = []
output_datetime = []
for t in output_time:
  date = ref_date + datetime.timedelta(days=t)
  output_date.append(date.strftime('%Y %m %d %H %M'))
  output_datetime.append(date)
output_datetime = np.asarray(output_datetime,dtype='O')

# Create station list
station_list = []
for run in stations:
  for sta in stations[run]['name']:
    station_list.append(sta)
station_list = list(set(station_list))


nc_data = nc4.Dataset(bathy_file)
lon_data = nc_data.variables['lon'][:]
lat_data = nc_data.variables['lat'][:]
z_data = nc_data.variables['z'][:,:]
#plt.figure()
#plt.contourf(lon_data,lat_data,z_data)
#plt.colorbar()
#plt.savefig('bathy.png')

bathy = interpolate.RegularGridInterpolator((lat_data,lon_data),z_data)
station_depth = []
for sta in station_list:
  for run in stations:
    if sta in stations[run]['name']:
      ind = stations[run]['name'].index(sta)
      lon = stations[run]['lon'][ind]
      lat = stations[run]['lat'][ind]
      break
  pt = np.array([lat,lon])
  depth = bathy(pt)
  station_depth.append(depth)


station_list = [x for _,x in sorted(zip(station_depth,station_list))]
station_depth = sorted(station_depth)


# Get data from observation file at output times
obs_data = {}
station_list_data = []
station_depth_data = []
for sta in station_list:
  
  obs_file = obs_direc+sta+'_'+year+'.txt'
  if os.path.isfile(obs_file):
    obs_data[sta] = plot_points.read_station_data(obs_file,output_date,variables)

    mask = np.isnan(obs_data[sta][variable])
    idx, = np.where(mask == False)

    if len(idx) > 1:
      station_list_data.append(sta)
      ind = station_list.index(sta)
      station_depth_data.append(station_depth[ind][0])


nsta = len(station_list_data)
for i in range(nsta):
  print station_list_data[i],station_depth_data[i]


xv = np.arange(nsta)
scatters = []
labels = []	

fig = plt.figure(figsize=[24,4])
ax = fig.add_subplot(111)
ax2 = ax.twinx()
dp, = ax2.plot(xv,station_depth_data,color='silver')

for run in direc:

  yv = np.zeros(nsta)
  yv.fill(np.nan)
  for sta in station_list_data:
    print sta
  
    if sta in stations[run]['name']:
    
      ind = stations[run]['name'].index(sta)
      lon = stations[run]['lon'][ind]
      lat = stations[run]['lat'][ind]
   
      # Only plot stations with modeled data (wet points)
      if len(np.unique(data[run][variable][:,ind])) > 1:
   
        mask = np.isnan(obs_data[sta][variable])
        idx, = np.where(mask == False)
        x = data[run][variable][:,ind][idx]
        y = obs_data[sta][variable][idx]
        
        # Only plot stations with observed data
        if len(idx) > 1:
   
          # Calculate regression coeffieint
#          slope,b,r,p,err = stats.linregress(x,y)
#          r2 = r**2
#          print r2

          r2 = np.mean(np.absolute(np.divide(y-x,y)))*100.0
   
          ind = station_list_data.index(sta)
          yv[ind] = r2
    
  
  sc = ax.scatter(xv,yv,marker='o',zorder=10)
  scatters.append(sc)
  labels.append(run)

scatters.append(dp)
labels.append('Depth at station')
lgd = plt.legend(scatters,labels,loc=2)
ax.set_xticks(xv)
ax.set_xticklabels(station_list_data,rotation='vertical')      
ax.set_xlabel('Station ID')
#ax.set_ylabel('r-squared')
ax.set_ylabel('Mean absolute percentage error')
ax2.set_ylabel('Depth (m)')
fig.savefig('rsquared.png',bbox_inches='tight',dpi=400)    
