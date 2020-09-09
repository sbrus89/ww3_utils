import plot_points
from scipy import stats
from scipy import interpolate
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import glob
import datetime
import numpy as np
import os
import netCDF4 as nc4
from collections import OrderedDict
import pandas as pd
import joypy
import math
import pickle

direc = OrderedDict()
#direc['1/2_degree'] = '/users/sbrus/scratch4/WW3_testing/glo_30m/post_processing/june/model_data/points/'
#direc['2_degree'] = '/users/sbrus/scratch4/WW3_testing/glo_2d/post_processing/june/model_data/points/'
#direc['unstr 4000 UOST'] = '/users/sbrus/scratch4/WW3_testing/unst_depth_4000_UOST/post_processing/june/model_data/points/'

direc['2 degree'] = '/users/sbrus/scratch4/WW3_testing/glo_2d/post_processing/june-october/model_data/points/'
direc['unstr 4000 UOST'] = '/users/sbrus/scratch4/WW3_testing/unst_depth_4000_UOST/post_processing/june-october/model_data/points/'
direc['1/2 degree'] = '/users/sbrus/scratch4/WW3_testing/glo_30m/post_processing/june-october/model_data/points/'

obs_direc = './obs_data/'
variables = plot_points.variables
variable = 'hs'
#variable = 'fp'
#variable = 'th1p'
year = '2005'
bathy_file = '/usr/projects/regionalclimate/COMMON_MPAS/ocean/grids/bathymetry_database/SRTM15_plus/earth_relief_15s.nc'

stations_exclude = ['44035','46228','44039','44013','46026','51202','44033','46216','46087','46081','46053','46054']

# Read in modeled data
data = {}
stations = {}
for run in direc:
  print(run)
  print('--------------------------------------------------')
  wav_files = [sorted(glob.glob(direc[run]+'ww3*_tab.nc'))]
  data[run], stations[run] = plot_points.read_point_files(wav_files,variables)

# Create station list
station_list = []
for run in stations:
  for sta in stations[run]['name']:
    station_list.append(sta)
station_list = list(set(station_list))
station_list_excluded = []
for sta in station_list:
  if sta in stations_exclude:
    continue
  else: 
    station_list_excluded.append(sta)
station_list = station_list_excluded


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
  print('interpolating '+sta)
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

# Find overall date range of data
date_min = '3000 01 01 00 00'
date_max = '1000 01 01 00 00'
frmt = '%Y %m %d %H %M'
for run in data:
  if data[run]['datetime'][0]  < datetime.datetime.strptime(date_min,frmt):
    date_min = data[run]['date'][0]
  if data[run]['datetime'][-1] > datetime.datetime.strptime(date_max,frmt):
    date_max = data[run]['date'][-1]
  runname = run
print(date_min,date_max)


# Get data from observation file at output times
obs_data = {}
station_list_data = []
station_depth_data = []
for sta in station_list:
  print('reading '+sta)
 
  obs_file = obs_direc+sta+'_'+year+'.txt'
  if os.path.isfile(obs_file):
    obs_data[sta] = plot_points.read_station_data(obs_file,date_min,date_max,variables,data[runname]['datetime'])

    mask = np.isnan(obs_data[sta][variable])
    idx, = np.where(mask == False)

    if len(idx) > 1:
      station_list_data.append(sta)
      ind = station_list.index(sta)
      station_depth_data.append(station_depth[ind][0])


nsta = len(station_list_data)
for i in range(nsta):
  print(station_list_data[i],station_depth_data[i])


xv = np.arange(nsta)
scatters = []
labels = []	

metrics = ['r2','mape','smape','rmse','mean error','abs error','rel bias','skill']
nmetrics = len(metrics)

fig = plt.figure(figsize=[24,5*nmetrics])
ax = []
for i in range(nmetrics):
  ax.append(fig.add_subplot(nmetrics,1,i+1))
  ax2 = ax[-1].twinx()
  dp, = ax2.plot(xv,station_depth_data,color='silver')
  ax2.set_ylabel('Depth (m)')

errors = {}
errors['station'] = []
errors_label = []
errors_abs = {}
errors_abs['station'] = []
errors_abs_label = []
errors_rel = {}
errors_rel['station'] = []
errors_rel_label = []
errors_absrel = {}
errors_absrel['station'] = []
errors_absrel_label = []

for irun,run in enumerate(direc):
  
  errors[run] = []
  errors_abs[run] = []
  errors_rel[run] = []
  errors_absrel[run] = []

  yv = np.zeros((nsta,nmetrics))
  yv.fill(np.nan)
  for ista,sta in enumerate(station_list_data):
    print(sta)
  
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
        
        if variable == 'fp':
          x = 1.0/x 
        if variable == 'th1p':
          d = np.degrees(np.arctan2(np.sin(np.radians(x-y)),np.sin(np.radians(x-y))))
        else:
          d = x-y


        
        # Only plot stations with observed data
        if len(idx) > 1:

          ind = station_list_data.index(sta)
          for i,metric in enumerate(metrics):
   
            if metric == 'r2' and variable != 'th1p':
              slope,b,r,p,err = stats.linregress(x,y)
              e = r**2
              yv[ind,i] = e 

            if metric == 'mape' and variable != 'th1p':
              e = np.mean(np.absolute(np.divide(d,y)))*100.0
              yv[ind,i] = e 

            if metric == 'smape' and variable != 'th1p':
              e = np.mean(np.divide(np.absolute(d),np.absolute(y)+np.absolute(x)))*100.0
              yv[ind,i] = e 

            if metric == 'rmse':
              e = np.sqrt(np.mean(np.square(d)))
              yv[ind,i] = e 

            if metric == 'mean error':
              e = np.mean(d)
              yv[ind,i] = e 

            if metric == 'abs error':
              e = np.mean(np.absolute(d))
              yv[ind,i] = e 

            if metric == 'rel bias' and variable != 'th1p':
              e = np.mean(d)/np.mean(y)
              yv[ind,i] = e 

            if metric == 'skill' and variable != 'th1p':
              e = 1.0 - np.sum(np.square(d))/np.sum(np.square(np.absolute(x-np.mean(y))+np.absolute(y-np.mean(y))))
              yv[ind,i] = e 

          for x in np.nditer(d):
            errors[run].append(x)
            if irun == 0:
              errors['station'].append(ista)
              if sta not in errors_label:
                errors_label.append(str(sta))
          for x in np.nditer(np.absolute(d)):
            errors_abs[run].append(x)
            if irun == 0:
              errors_abs['station'].append(ista)
              if sta not in errors_abs_label:
                errors_abs_label.append(str(sta))
          for x in np.nditer(100.0*np.divide(d,y)):
            if np.isnan(x) or np.isinf(x):
              continue
            else:
              errors_rel[run].append(x)
              if irun == 0:
                errors_rel['station'].append(ista)
                if sta not in errors_rel_label:
                  errors_rel_label.append(str(sta))
          for x in np.nditer(100.0*np.absolute(np.divide(d,y))):
            if np.isnan(x) or np.isinf(x):
              continue
            else:
              errors_absrel[run].append(x)
              if irun == 0:
                errors_absrel['station'].append(ista)
                if sta not in errors_absrel_label:
                  errors_absrel_label.append(str(sta))
              
                 
  for i in range(nmetrics):
    sc = ax[i].scatter(xv,yv[:,i],marker='o',zorder=10)
    if i == 0:
      scatters.append(sc)
      labels.append(run)

  with open(run.replace(' ','_').replace('/','-')+'.pickle','wb') as f:
    obj = [run,xv,yv,metrics]
    pickle.dump(obj,f,protocol=pickle.HIGHEST_PROTOCOL)

scatters.append(dp)
labels.append('Depth at station')
lgd = ax[0].legend(scatters,labels,loc=2)
for i,metric in enumerate(metrics):
  ax[i].set_xticks(xv)
  ax[i].set_xticklabels(station_list_data,rotation='vertical')      
  ax[i].set_xlabel('Station ID')
  if metric == 'r2':
    ax[i].set_ylabel('r-squared')
  if metric == 'mape':
    ax[i].set_ylabel('Mean absolute percentage error')
  if metric == 'smape':
    ax[i].set_ylabel('symmetric mean absolute percentage error')
  if metric == 'rmse':
    ax[i].set_ylabel('root mean square error')
  if metric == 'mean error':
    ax[i].set_ylabel('mean error')
  if metric == 'abs error':
    ax[i].set_ylabel('absolute error')
  if metric == 'rel bias':
    ax[i].set_ylabel('relative bias')
  if metric == 'skill':
    ax[i].set_ylabel('model skill')
fig.savefig('rsquared.png',bbox_inches='tight',dpi=400)    

for key in errors:
  if key != 'station' and key != 'label':
    errors[key] = np.asarray(errors[key],dtype=np.float64)  
    errors_abs[key] = np.asarray(errors_abs[key],dtype=np.float64)  
    errors_rel[key] = np.asarray(errors_rel[key],dtype=np.float64)  
    errors_absrel[key] = np.asarray(errors_absrel[key],dtype=np.float64)  


with open('errors.pickle','wb') as f:
  errors_obj = [errors,errors_label]
  pickle.dump(errors_obj,f,protocol=pickle.HIGHEST_PROTOCOL)
with open('errors_abs.pickle','wb') as f:
  errors_abs_obj = [errors_abs,errors_abs_label]
  pickle.dump(errors_abs_obj,f,protocol=pickle.HIGHEST_PROTOCOL)
with open('errors_rel.pickle','wb') as f:
  errors_rel_obj = [errors_rel,errors_rel_label]
  pickle.dump(errors_rel_obj,f,protocol=pickle.HIGHEST_PROTOCOL)
with open('errors_absrel.pickle','wb') as f:
  errors_absrel_obj = [errors_absrel,errors_absrel_label]
  pickle.dump(errors_absrel_obj,f,protocol=pickle.HIGHEST_PROTOCOL)

columns = list(errors.keys()) 
columns.remove('station')
fsize = (4,18)
opac = 0.5

df = pd.DataFrame.from_dict(errors)
plt.figure(figsize=fsize)
fig, axes = joypy.joyplot(df, column=columns,by='station',kind='normalized_counts',bins=30,alpha=opac, figsize=fsize, x_range=[-2,2],labels=errors_label, legend=True)
plt.title('Errors')
plt.savefig('joyplot_errors.png',bbox_inches='tight',dpi=400)    

df = pd.DataFrame.from_dict(errors_abs)
plt.figure(figsize=fsize)
fig, axes = joypy.joyplot(df, column=columns,by='station',kind='normalized_counts',bins=30,alpha=opac, figsize=fsize, x_range=[0,2], labels=errors_abs_label, legend=True)
plt.title('Absolute errors')
plt.savefig('joyplot_errors_abs.png',bbox_inches='tight',dpi=400)    

df = pd.DataFrame.from_dict(errors_rel)
plt.figure(figsize=fsize)
fig, axes = joypy.joyplot(df, column=columns,by='station',kind='normalized_counts',bins=30,alpha=opac, figsize=fsize, x_range=[-200,200],labels=errors_rel_label, legend=True)
plt.title('Relative errors')
plt.savefig('joyplot_errors_rel.png',bbox_inches='tight',dpi=400)    

df = pd.DataFrame.from_dict(errors_absrel)
plt.figure(figsize=fsize)
fig, axes = joypy.joyplot(df, column=columns,by='station',kind='normalized_counts',bins=30,alpha=opac, figsize=fsize, x_range=[0,200],labels=errors_absrel_label, legend=True)
plt.title('Absolute relative errors')
plt.savefig('joyplot_errors_absrel.png',bbox_inches='tight',dpi=400)    
