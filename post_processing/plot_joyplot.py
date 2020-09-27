import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib.cm as cm
import matplotlib.colors as colors
import numpy as np
import pandas as pd
import joypy
import pickle
import plot_points
import pprint
import shapely.geometry
from geometric_features import read_feature_collection
import netCDF4 as nc4 
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from adjustText import adjust_text
import subprocess

rc('font',**{'family':'serif','serif':['Times New Roman']})

average_error_files = ['2_degree_UOST_PR1.pickle','unstr_4000_UOSTcor_edit2.pickle','1-2_degree_UOST_PR1.pickle']

stations = plot_points.read_station_file('stations.txt')

filename = 'etopo1.nc' 
nc_data = nc4.Dataset(filename, "r")
lon_data = nc_data.variables['lon'][:]
lat_data = nc_data.variables['lat'][:]
bathy_data = nc_data.variables['z'][:]*-1.0
idx = np.where(lon_data > 180.0)
lon_data[idx] = lon_data[idx] - 360.0
idx = np.argsort(lon_data)
lon_data = lon_data[idx]
bathy_data = bathy_data[:,idx]

with open('errors.pickle','rb') as f:
  errors,errors_label = pickle.load(f)
with open('errors_abs.pickle','rb') as f:
  errors_abs,errors_abs_label = pickle.load(f)
with open('errors_rel.pickle','rb') as f:
  errors_rel,errors_rel_label = pickle.load(f)
with open('errors_absrel.pickle','rb') as f:
  errors_absrel,errors_absrel_label = pickle.load(f)

average_errors = []
for fname in average_error_files:
  with open(fname,'rb') as f:
    run,xv,yv,metrics,station_list_data,station_depth_data = pickle.load(f)
    average_errors.append({'run_name':run,'xv':xv,'yv':yv,'metrics':metrics,'stations':station_list_data,'depth':station_depth_data})

columns = list(errors.keys()) 
columns.remove('station')
fsize = (3,9)
opac = 0.5

pprint.pprint(errors_label)
pprint.pprint(len(errors_label))
pprint.pprint(max(errors['station']))

fc = read_feature_collection('NDBC_regions.geojson')

regions = [ 
           'East_Coast_South',
           'East_Coast_North',
           'HI',
           'AK',
           'Gulf',
           'West_Coast_South',
           'West_Coast_North',
           'Caribbean_Atlantic'
         ]

for region in regions:
  for feature in fc.features:


    ###########################################
    #  Find stations in region
    ###########################################  
    if feature['properties']['name'] == region:
      pass
    else:
      continue
    print(feature['geometry'])
  
    correct_dateline = False
    if np.min(feature['geometry']['coordinates'][0]) < -180.0:
      correct_dateline = True
      print('corrected dateline')
    shape = shapely.geometry.Polygon(feature['geometry']['coordinates'][0])
  
    stations_include = {}
    stations_include['idx'] = []
    stations_include['name'] = []
    stations_include['coords'] = [] 
    for i,sta in enumerate(errors_label):
      idx = stations['name'].index(sta)
  
      if correct_dateline and stations['lon'][idx] > 0.0:
        lon = stations['lon'][idx] - 360.0
      else: 
        lon = stations['lon'][idx]
  
      point = shapely.geometry.Point(lon,stations['lat'][idx])
      test = shape.contains(point)
      if test == True:
        print(sta,feature['properties']['name'])
        stations_include['idx'].append(i)
        stations_include['name'].append(sta)
        if correct_dateline and stations['lon'][idx] < 0.0:
          lon = stations['lon'][idx] + 360.0
        else:
          lon = stations['lon'][idx]
        stations_include['coords'].append([lon,stations['lat'][idx]])
  
    pprint.pprint(stations_include)
    print(len(stations_include['name']))
   
  ###########################################
  # Joy plot 
  ###########################################  
  print('creating joy plot') 
  
  df = pd.DataFrame.from_dict(errors)
  df_new = df[df['station'].isin(stations_include['idx'])]
  
  plt.figure(figsize=fsize)
  if region == 'West_Coast_North':
    overlap = 0.5
  elif region == 'West_Coast_South':
    overlap = 0.8
  elif region == 'Gulf':
    overlap = 0.5
  elif region == 'AK':
    overlap = 0.7
  #elif region == 'HI':
  #  overlap = 0.5
  elif region == 'East_Coast_South':
    overlap = 0.8
  elif region == 'East_Coast_North':
    overlap = 0.5
  else:
    overlap = 1.0
  fig, axes = joypy.joyplot(df_new, grid=True, column=columns,by='station',kind='normalized_counts',bins=30,alpha=opac, figsize=fsize, x_range=[-2,2],labels=stations_include['name'], legend=False, xlabelsize=12,ylabelsize=12,overlap=overlap)
  axes[-1].set_title('(a)',loc='center',fontsize=18)
  axes[-1].set_xlabel('Error (m)',fontsize=12)
  axes[-1].yaxis.set_visible(True)
  axes[-1].yaxis.set_ticks([])
  axes[-1].set_ylabel('Station ID',fontsize=12,labelpad=55)
  plt.savefig('joyplot_errors_'+region+'.png',bbox_inches='tight',dpi=400)    
  print('  done')

  ###########################################
  # Station locations 
  ###########################################  
  
  print('creating station location plot')
  stations_include['coords'] = np.array(stations_include['coords'])
  crs0 = ccrs.PlateCarree(central_longitude=0)
  crs180 = ccrs.PlateCarree(central_longitude=180)
  if region == 'West_Coast_South':
    extent_pad = 0.5
  elif region == 'West_Coast_North':
    extent_pad = 2.0
  elif region == 'East_Coast_North':
    extent_pad = 1.0
  elif region == 'AK':
    extent_pad = 4.0
  elif region == 'Gulf':
    extent_pad = 2.0
  elif region == 'HI':
    extent_pad = 1.0
  else:
    extent_pad = 5.0
  lon_min = np.min(stations_include['coords'][:,0])-extent_pad 
  lon_max = np.max(stations_include['coords'][:,0])+extent_pad
  lat_min = np.min(stations_include['coords'][:,1])-extent_pad
  lat_max = np.max(stations_include['coords'][:,1])+extent_pad
  lon_mid = 0.5*(lon_min+lon_max)
  lat_mid = 0.5*(lat_min+lat_max)
  pm = 0.5*max(lon_max-lon_min,lat_max-lat_min)
  lon_min = lon_mid-pm
  lon_max = lon_mid+pm
  lat_min = lat_mid-pm
  lat_max = lat_mid+pm
  extent = [lon_min,lon_max,lat_min,lat_max]
  fig = plt.figure(figsize=[6,6])
  ax = fig.add_subplot(1,1,1,projection=crs180)
  lon_idx, = np.where((lon_data >= lon_min) & (lon_data <= lon_max))
  lat_idx, = np.where((lat_data >= lat_min) & (lat_data <= lat_max))
  lon_region = lon_data[lon_idx]
  lat_region = lat_data[lat_idx]
  latlon_idx = np.ix_(lat_idx, lon_idx)
  bathy_region = bathy_data[latlon_idx]
  #ds = 50                                # Downsample
  #dsx = np.arange(0, lon_data.size, ds)  # bathy data
  #dsy = np.arange(0, lat_data.size, ds)  # to speed up
  #dsxy = np.ix_(dsy, dsx)                # plotting
  colormap = cm.get_cmap('Blues')
  minval = 0.25
  maxval = 1.0 
  cmap = colors.LinearSegmentedColormap.from_list(
       'trunc({n},{a:.2f},{b:.2f})'.format(n=colormap.name, a=minval, b=maxval),
       colormap(np.linspace(minval, maxval, 100)))

  levels = np.linspace(0,8000,50)
  trans = crs0
  ax.contourf(lon_region,lat_region,bathy_region,levels=levels,transform=trans,cmap=cmap)
  ax.set_extent(extent)
  ax.set_xticks(np.round_(np.linspace(lon_min,lon_max,5),1),crs=ccrs.PlateCarree())
  ax.set_yticks(np.round_(np.linspace(lat_min,lat_max,5),1),crs=ccrs.PlateCarree())
  lon_formatter = LongitudeFormatter(zero_direction_label=True)
  lat_formatter = LatitudeFormatter()
  ax.xaxis.set_major_formatter(lon_formatter)
  ax.yaxis.set_major_formatter(lat_formatter)
  ax.set_xlabel('Longitude',fontsize=12)
  ax.set_ylabel('Latitiude',fontsize=12)
  annotations = []
  for i,sta in enumerate(stations_include['name']):
    at_x, at_y = crs180.transform_point(stations_include['coords'][i,0],stations_include['coords'][i,1],crs0)
    annotations.append(ax.text(at_x,at_y,sta,fontsize=14,fontweight='bold'))
  adjust_text(annotations,  arrowprops=dict(arrowstyle="->", color='k', lw=0.7,alpha=0.5),precision=0.001)
  
  for i,sta in enumerate(stations_include['name']):
    x,y = crs180.transform_point(stations_include['coords'][i,0],stations_include['coords'][i,1],crs0)
    ax.scatter(x,y,transform=crs180,c='C3')
  ax.add_feature(cfeature.LAND)
  ax.add_feature(cfeature.COASTLINE,edgecolor='gray')
  ax.set_title('(c)',loc='center',fontsize=18)
  
  plt.savefig('station_locations_'+region+'.png',bbox_inches='tight',dpi=400)
  print('  done')
  
  ###########################################
  # Average error plot
  ###########################################  

  print('creating average error plot')
  selected_metrics = ['rmse']
  nmetric = len(average_errors[0]['metrics'])
  nrun = len(average_errors)
  fig = plt.figure(figsize=(6,3*len(selected_metrics)))
  idx = np.array(stations_include['idx'],dtype=np.int32)
  xlabels = []
  for sta in stations_include['idx']:
    xlabels.append(average_errors[0]['stations'][sta])
  station_depths = np.array(average_errors[0]['depth'])
  nsta = idx.size
  xv = np.arange(nsta)
  n = 0
  scatters = []
  labels = []
  for i in range(nmetric):
    metric = average_errors[0]['metrics'][i]
    if metric in selected_metrics:
      n = n + 1
    else:
      continue
    ax = fig.add_subplot(len(selected_metrics),1,n)
  
    for j in range(nrun):
      sc = ax.scatter(xv,average_errors[j]['yv'][idx,i],marker='o',zorder=10,alpha=0.6)
      if n == 1:
        scatters.append(sc)
        labels.append(average_errors[j]['run_name'])
    ax.set_xticks(xv)
    ax.set_xticklabels(xlabels,rotation='vertical')
    ax.tick_params(axis='both',which='major',labelsize=12)
    ax2 = ax.twinx()
    dp, = ax2.plot(xv,station_depths[idx],color='silver')
    ax2.set_ylabel('Depth (m)',fontsize=12)
    
    if metric == 'r2':
      ylabel = 'r-squared'
    if metric == 'mape':
      ylabel = 'Mean absolute percentage error'
    if metric == 'smape':
      ylabel = 'symmetric mean absolute percentage error'
    if metric == 'rmse':
      ylabel = 'RMSE (m)'
    if metric == 'mean error':
      ylabel = 'mean error (m)'
    if metric == 'abs error':
      ylabel = 'absolute error (m)'
    if metric == 'rel bias':
      ylabel = 'relative bias (m)'
    if metric == 'skill':
      ylabel='model skill'
    ax.set_ylabel(ylabel,fontsize=12)
    ax.set_xlabel('Station ID',fontsize=12)
  ax.set_title('(b)',loc='center',fontsize=18)
  plt.savefig('average_errors_'+region+'.png',bbox_inches='tight',dpi=400)
  print('  done')
  
  ###########################################
  # Legend 
  ###########################################  
  
  print('creating legend')
  fig = plt.figure(figsize=(6,1))
  for i in range(len(labels)):
    if labels[i].find('1/2 degree') >= 0:
      labels[i] = '1/2 degree structured'
    elif labels[i].find('2 degree') >= 0:
      labels[i] = '2 degree structured'
    elif labels[i].find('unst') >= 0:
      labels[i] = 'unstructured'
  fig.legend(scatters,labels,loc='center',fontsize=12,ncol=3,fancybox=True)
  plt.savefig('legend.png',bbox_inches='tight',dpi=400)    
  print('  done')
  
  ###########################################
  # Concatenate figures 
  ###########################################  
  
  print('concatenating figures')
  subprocess.call('convert -append average_errors_'+region+'.png station_locations_'+region+'.png legend.png out.png',shell=True)
  subprocess.call('convert +append joyplot_errors_'+region+'.png out.png '+region+'.png',shell=True)
  print('  done')


#df = pd.DataFrame.from_dict(errors_abs)
#plt.figure(figsize=fsize)
#fig, axes = joypy.joyplot(df, column=columns,by='station',kind='normalized_counts',bins=30,alpha=opac, figsize=fsize, x_range=[0,2], labels=errors_abs_label, legend=False)
#plt.title('Absolute errors')
#plt.savefig('joyplot_errors_abs.png',bbox_inches='tight',dpi=400)    
#
#df = pd.DataFrame.from_dict(errors_rel)
#plt.figure(figsize=fsize)
#fig, axes = joypy.joyplot(df, column=columns,by='station',kind='normalized_counts',bins=30,alpha=opac, figsize=fsize, x_range=[-200,200],labels=errors_rel_label, legend=False)
#plt.title('Relative errors')
#plt.savefig('joyplot_errors_rel.png',bbox_inches='tight',dpi=400)    
#
#df = pd.DataFrame.from_dict(errors_absrel)
#plt.figure(figsize=fsize)
#fig, axes = joypy.joyplot(df, column=columns,by='station',kind='normalized_counts',bins=30,alpha=opac, figsize=fsize, x_range=[0,200],labels=errors_absrel_label, legend=False)
#plt.title('Absolute relative errors')
#plt.savefig('joyplot_errors_absrel.png',bbox_inches='tight',dpi=400)    
