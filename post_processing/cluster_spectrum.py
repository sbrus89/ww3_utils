import glob
import yaml
import os
import pprint
import time
import socket
import sys
import numpy as np
import netCDF4
from matplotlib import cm,ticker
from matplotlib.transforms import Bbox
from matplotlib.patches import Rectangle
from matplotlib.colors import ListedColormap
from matplotlib.colors import LogNorm 
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy import interpolate
from scipy import sparse
from scipy import stats
import directional_spectrum
from sklearn.cluster import KMeans
from sklearn.cluster import DBSCAN 
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import SpectralClustering
from sklearn.mixture import GaussianMixture
from sklearn.decomposition import PCA
from sklearn.metrics import r2_score
np.set_printoptions(threshold=sys.maxsize,linewidth=1000)

# Don't use interactive plotting on IC machines
host = socket.gethostname()
if host[0:2] == 'gr' or host[0:2] == 'ba':
  plt.switch_backend('agg')
  interactive = False
else:
  interactive = True

# Read in config file 
pwd = os.getcwd()
config_file = pwd+'/cluster_spectrum.config'
if os.path.exists(config_file):
    # Read in config file
    f = open(config_file)
    cfg = yaml.load(f,Loader=yaml.Loader)
    pprint.pprint(cfg)
    run_direcs = cfg['run_direcs']
    month = cfg['month']
else:
  run_direcs = ['./test_model_data/']
  month = 7

##########################################################################
##########################################################################
# User-defined options 
##########################################################################
##########################################################################

k_vals = [16]
#k_vals = [16,25,36] 
cluster_methods = ['kmeans']
#cluster_methods = ['agglomerative','kmeans','mixture']
normalize = True
error_metrics = ['RMSE','L2','SMAPE','R2'] 
lat_min = -65.0
lat_max = 65.0
log_scale = False
#interactive = False

##########################################################################
##########################################################################
# Input functions
##########################################################################
##########################################################################

def get_averaged_model_spectra(run_direcs,month,normalize):

  # Read in and average model results
  print('Averaging station spectra')
  run_files = []
  for direc in run_direcs:
    run_files.extend(sorted(glob.glob(direc+'spec.????'+str(month).zfill(2)+'*T??Z_spec.nc')))
  if len(run_files) > 0:
    run_data,run_stations = directional_spectrum.read_model_data(run_files,'average')

  # Find only grid stations only (no NDBC bouys)
  print('Eliminating NDBC and ice covered stations')
  m = 0
  for sta in run_stations['name']:
    if sta.find(',') < 0:
      m = m+1
  run_stations['name'] = run_stations['name'][m:]
  run_stations['lon'] = run_stations['lon'][m:]
  run_stations['lat'] = run_stations['lat'][m:]
  nstations = len(run_stations['name'])

  # Get spectral data for grid stations
  one,total_stations,nfrequency,ndirection = run_data['spectrum'].shape
  spectrum = np.reshape(run_data['spectrum'][0,m:,:,:],(nstations,nfrequency*ndirection))
  print(spectrum.shape)
  
  # Get rid of stations covered with ice
  nonice_stations = []
  ice_stations = []
  for sta in range(nstations):
    row_max = np.amax(spectrum[sta,:])
    if row_max > 0.:
      nonice_stations.append(sta)
    else:
      ice_stations.append(sta)
  spectrum = np.delete(spectrum,ice_stations,0)
  print(spectrum.shape)

  run_stations['name'] = [run_stations['name'][i] for i in nonice_stations]
  run_stations['lon'] = run_stations['lon'][nonice_stations]
  run_stations['lat'] = run_stations['lat'][nonice_stations]
  nstations = len(run_stations['name'])

  # Get rid of stations in certain locations 
  in_stations = []
  out_stations = []
  for sta,lat in enumerate(run_stations['lat']):
    lon = run_stations['lon'][sta]
    if lat > lat_max or lat < lat_min:                                  # high latitudes
      out_stations.append(sta)
    elif lon > 0.0 and lon < 40.0 and lat > 30.0 and lat < 45:          # Mediterranean Sea
      out_stations.append(sta)
    elif lon > -100.0 and lon < -65.0 and lat > 50.0 and lat < lat_max: # Hudson Bay
      out_stations.append(sta)
    elif lon > 47.5 and lon < 56.5 and lat > 23.5 and lat < 30.0:       # Persian Gulf
      out_stations.append(sta)
    else: 
      in_stations.append(sta)
  spectrum = np.delete(spectrum,out_stations,0)
  print(spectrum.shape)

  run_stations['name'] = [run_stations['name'][i] for i in in_stations]
  run_stations['lon'] = run_stations['lon'][in_stations]
  run_stations['lat'] = run_stations['lat'][in_stations]
  nstations = len(run_stations['name'])
    

  # Normalize the spectra
  if normalize:
    print('Normalizing spectra')
    for sta in range(nstations):
      row_max = np.amax(spectrum[sta,:])
      row_min = np.amin(spectrum[sta,:])
      if row_max > 0.:
        spectrum[sta,:] = (spectrum[sta,:]-row_min)/(row_max-row_min)

  return run_data,run_stations,spectrum

##########################################################################
##########################################################################

def write_averaged_model_spectra(filename,run_data,run_stations,spectrum):
 
  nc_file = netCDF4.Dataset(filename,'w')

  nsta,nfreqdir = spectrum.shape
  one,total_stations,nfrequency,ndirection = run_data['spectrum'].shape
  nc_file.createDimension('nstations',nsta)
  nc_file.createDimension('ndirection',ndirection)
  nc_file.createDimension('nfrequency',nfrequency)
  nc_file.createDimension('nfreqdir',ndirection*nfrequency)

  freq_nc = nc_file.createVariable('freq','f8',('nfrequency',))
  direc_nc = nc_file.createVariable('direc','f8',('ndirection',)) 
  spec_nc = nc_file.createVariable('spec','f8',('nstations','nfreqdir'))
  lat_nc = nc_file.createVariable('lat','f8',('nstations',))
  lon_nc = nc_file.createVariable('lon','f8',('nstations',))
  
  freq_nc[:] = run_data['freq']
  direc_nc[:] = run_data['theta']
  spec_nc[:] = spectrum
  lon_nc[:] = run_stations['lon']
  lat_nc[:] = run_stations['lat']

  nc_file.close() 
   
##########################################################################
##########################################################################

def read_averaged_model_spectra(filename):
 
  nc_file = netCDF4.Dataset(filename,'r')

  run_data = {}
  run_stations = {}

  run_data['freq'] = nc_file.variables['freq'][:]
  run_data['theta'] = nc_file.variables['direc'][:]
  run_stations['lon'] = nc_file.variables['lon'][:]
  run_stations['lat'] = nc_file.variables['lat'][:]
  spectrum = nc_file.variables['spec'][:]

  nc_file.close() 

  return run_data,run_stations,spectrum

##########################################################################
##########################################################################
# Clustering functions
##########################################################################
##########################################################################

def cluster_spectra(k,method,spectrum,run_stations):

  print('Performing '+method+' clustering')
  t0 = time.time()
  if method == 'kmeans':
    clustered = KMeans(n_clusters=k,random_state=0).fit(spectrum)
    spec = clustered.cluster_centers_
    labels = clustered.labels_
  elif method == 'pca': 
    clustered = PCA(n_components=k).fit(spectrum)
    spec = clustered.components_
  elif method == 'dbscan':
    clustered = DBSCAN(eps=4.0,min_samples=1).fit(spectrum)
    labels = clustered.labels_
    k = np.amax(labels)+1
    idx, = np.where(labels==-1)
    print(np.amax(labels))
    print(idx.shape)
    spec = np.zeros((k,spectrum.shape[1]))
    print(spec.shape)
    print(spectrum.shape)
    for i in range(k):
      idx, = np.where(labels==i)
      spec[i,:] = np.mean(spectrum[idx,:],axis=0)
  elif method == 'mixture':
    clustered = GaussianMixture(n_components=k)
    spec = clustered.fit(spectrum).means_
    labels = clustered.predict(spectrum)
  elif method == 'agglomerative':
    #station_connectivity = compute_connectivity(run_stations['name'])
    #clustered = AgglomerativeClustering(n_clusters=k,linkage='ward',connectivity=station_connectivity).fit(spectrum)
    clustered = AgglomerativeClustering(n_clusters=k,linkage='ward').fit(spectrum)
    labels = clustered.labels_
    spec = calculate_cluster_averages(k,labels,spectrum) 
  elif method == 'spectral':
    #clustered = SpectralClustering(n_clusters=k,affinity=calculate_R2,random_state=0).fit(spectrum)
    clustered = SpectralClustering(n_clusters=k,affinity='linear',random_state=0).fit(spectrum)
    labels = clustered.labels_
    spec = calculate_cluster_averages(k,labels,spectrum) 
  t1 = time.time()
  print("elapsed time = ",t1-t0)
    
  return labels,spec,k

##########################################################################
##########################################################################

def compute_connectivity(stations):

  ddeg = 3.0
  
  lon_min = -180.0
  lon_max =  180.0
  lon_vec = np.arange(lon_min,lon_max,ddeg)
  nlon = lon_vec.size
  
  lat_min = -80.0
  lat_max =  80.0
  lat_vec = np.arange(lat_min,lat_max,ddeg)
  nlat = lat_vec.size
  
  lon_grid,lat_grid = np.meshgrid(lon_vec,lat_vec)
  lon_pts = np.ravel(lon_grid)
  lat_pts = np.ravel(lat_grid)
  station_grid = np.zeros(lon_grid.shape).astype(np.int32)-1    
  nstations = len(stations)
  station_connectivity = np.zeros((nstations,nstations)).astype(np.int32)
                                                                                                                                    
  n = 0
  for sta in stations:
    if sta.find(',') > 0:
      lon = float(sta.split(',')[0])
      lat = float(sta.split(',')[1])
      i = np.where(lat_vec == lat)
      j = np.where(lon_vec == lon)
      station_grid[i,j] = n
      n = n+1

  for i in range(0,nlat):
    for j in range(0,nlon):
      indicies= [[i-1,j],[i+1,j],[i,j-1],[i,j+1],[i-1,j+1],[i-1,j-1],[i+1,j+1],[i+1,j-1]]
      neighbors = []
      for idx in indicies:
        if idx[1] == nlon: 
          idxj = 0
        else:
          idxj = idx[1]
        if idx[0] < nlat-1 and idx[0] > 0:
          idxi = idx[0]
          if station_grid[idxi,idxj] > 0 :
            neighbors.append(station_grid[idxi,idxj])
      n = station_grid[i,j]
      for m in neighbors:
        station_connectivity[n,m] = 1
  station_connectivity = sparse.csr_matrix(station_connectivity)
  return station_connectivity
  
##########################################################################
##########################################################################

def calculate_cluster_averages(k,labels,model_spectra):

  cluster_averages = np.zeros((k,model_spectra.shape[1]))
  print(cluster_averages.shape)
  print(model_spectra.shape)
  for i in range(k):
    idx, = np.where(labels==i)
    cluster_averages[i,:] = np.mean(model_spectra[idx,:],axis=0)

  return cluster_averages

##########################################################################
##########################################################################
# Cluster plotting functions
##########################################################################
##########################################################################

def plot_station_scatter(k,method,labels,station_lon,station_lat,theta,freq,clustered_spectra,model_spectra,cmap):

  if method == 'pca':
    return

  print('Plotting clustered station labels')

  # Show actual and clustered spectra for selected point
  def onpick(event):
    ind = event.ind[0]
    print(event.ind)
  
    i = labels[ind]
    ndir = 37
    nfreq = 48

    nfrequency = freq.size
    ndirection = theta.size
     
    # Plot cluster center    
    spectrum_cluster = np.reshape(clustered_spectra[i,:],(nfrequency,ndirection)).T  
    Theta,Freq,spec_interp = directional_spectrum.interp_model_spectrum(theta,freq,spectrum_cluster,ndir,nfreq) 
    fig = plt.figure()
    ax = fig.add_subplot(1,2,1,polar=True)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    cax = ax.contourf(Theta,Freq,spec_interp,30)
    cb = fig.colorbar(cax)
    fig.tight_layout()
    rect = [Rectangle([0.0,0.0],0.5,1,facecolor=cmap(i),edgecolor='none',zorder=-1,transform=fig.transFigure)]
    fig.patches.extend(rect)
    ax.set_title('Custer Center')

    # Plot modelled station spectrum
    spectrum_model = np.reshape(model_spectra[ind,:],(nfrequency,ndirection)).T  
    Theta,Freq,spec_interp = directional_spectrum.interp_model_spectrum(theta,freq,spectrum_model,ndir,nfreq) 
    ax = fig.add_subplot(1,2,2,polar=True)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    cax = ax.contourf(Theta,Freq,spec_interp,30)
    cb = fig.colorbar(cax)
    ax.set_title('Modeled Spectrum')
    plt.plot()
    plt.show()

  # Plot scatterplot of clustered stations
  fig = plt.figure(figsize=[16.0,6.0])
  central_longitude = 205.0
  crs = ccrs.PlateCarree(central_longitude)
  ax = fig.add_subplot(111,projection=crs)
  ax.set_extent([0.0,359.99,np.amin(station_lat),np.max(station_lat)],crs=crs)
  ax.add_feature(cfeature.LAND)
  ax.add_feature(cfeature.LAKES)
  ax.add_feature(cfeature.COASTLINE)
  plt.title(method+' clusters for k='+str(k))
  if interactive:
    cax = ax.scatter(station_lon+360.0-central_longitude,station_lat,c=labels,cmap=cmap,picker=True,transform=crs)
    fig.canvas.mpl_connect('pick_event',onpick)
    cb = fig.colorbar(cax)
    fig.tight_layout()
    plt.plot()
    plt.show()
  else:
    cax = plt.scatter(np.mod(station_lon+360,360),station_lat,c=labels,cmap=cmap)
    cb = fig.colorbar(cax)
    fig.tight_layout()
  fig.savefig('clustered_spectrum_'+method+'_k='+str(k)+'.png',bbox_inches='tight')
  plt.close()

##########################################################################
##########################################################################

def plot_cluster_grid(k,method,clustered_spectra,theta,freq,cmap):
  print('Plotting grid of cluster centers')
  print(clustered_spectra.shape)
  
  ndir = 37
  nfreq = 48
  
  n = int(np.ceil(np.sqrt(k)))
  
  fig = plt.figure(figsize=[18,18])
  for i in range(k):
    spectrum_cluster = np.reshape(clustered_spectra[i,:],(freq.size,theta.size)).T  
    Theta,Freq,spec_interp = directional_spectrum.interp_model_spectrum(theta,freq,spectrum_cluster,ndir,nfreq) 
  
    # Plot spectrum
    ax = fig.add_subplot(n,n,i+1,polar=True)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    if log_scale:
      cax = ax.pcolormesh(Theta,Freq,spec_interp,norm=colors.LogNorm(vmin=1e-3,vmax=1.0))
      cb = fig.colorbar(cax,ticks=ticker.LogLocator(subs=range(10)),fraction=0.07, pad=0.1)
      ax.plot(np.linspace(0,2*np.pi,100), np.linspace(0,0.5,100), color='k', ls='none') 
      ax.grid()
    else:
      cax = ax.contourf(Theta,Freq,spec_interp,30)
      cb = fig.colorbar(cax,fraction=0.07, pad=0.1)
  
  # Plot color of cluster behind each subplot
  fig.tight_layout()
  h = 1.0/float(n)
  c = 0
  for j in range(n):
    for i in range(n):
     rect = [Rectangle([0.0+(i*h),1.0-(j*h)-h],h,h,facecolor=cmap(c),edgecolor='none',zorder=-1,transform=fig.transFigure)]
     c = c +1
     fig.patches.extend(rect)
  
  fig.savefig('k_spectra_'+method+'_k='+str(k)+'.png',bbox_inches='tight')
  plt.close()

##########################################################################
##########################################################################

def get_colormap(k):

  # Get colormap
  cmap1 = cm.get_cmap('tab20b')(range(20))
  cmap2 = cm.get_cmap('tab20c')(range(20))
  cmap = ListedColormap(np.vstack((cmap1,cmap2)))
  cmap = ListedColormap(cmap(range(k)))

  ## Plot colors
  #fig = plt.figure()
  #x = np.arange(k)
  #plt.scatter(x,x,c=x,cmap=cmap)
  #fig.savefig('colors.png',bbox_inches='tight')
  #plt.close()

  return cmap

##########################################################################
##########################################################################
# Error plotting functions
##########################################################################
##########################################################################

def plot_cluster_error(k,method,metric,labels,model_spectra,clustered_spectra,station_lon,station_lat):

  print('Plotting '+metric+' error between station spectra and cluster centers')
  nstations = model_spectra.shape[0]
  cluster_error = np.zeros(nstations)
  for sta in range(nstations):
    spectrum_model = model_spectra[sta,:]
    spectrum_cluster = clustered_spectra[labels[sta],:]
    #row_min = np.amin(spectrum_cluster)
    #row_max = np.amax(spectrum_cluster)
    #spectrum_cluster = (spectrum_cluster-row_min)/(row_max-row_min)
    if metric == 'RMSE':
      cluster_error[sta] = np.sqrt(np.mean(np.square(spectrum_model-spectrum_cluster)))
      vmax = 0.75*np.amax(cluster_error)
      vmin = np.amin(cluster_error)
    elif metric == 'L2':
      cluster_error[sta] = np.sqrt(np.sum(np.square(spectrum_model-spectrum_cluster)))
      vmax = 0.75*np.amax(cluster_error)
      vmin = np.amin(cluster_error)
    elif metric == 'SMAPE':
      cluster_error[sta] = np.mean(np.divide(np.absolute(spectrum_model-spectrum_cluster),np.absolute(spectrum_model)+np.absolute(spectrum_cluster)))
      vmax = 1.0
      vmin = 0.0
    elif metric == 'R2':
      slope,intercept,r_value,p_value,std_err = stats.linregress(spectrum_model,spectrum_cluster)
      cluster_error[sta] = r_value**2
      vmax = 1.0
      vmin = 0.0

  fig = plt.figure(figsize=[16.0,6.0])
  central_longitude = 205.0
  crs = ccrs.PlateCarree(central_longitude)
  ax = fig.add_subplot(111,projection=crs)
  ax.set_extent([0.0, 359.99, np.amin(station_lat), np.amax(station_lat)],crs=crs)
  ax.add_feature(cfeature.LAND)
  ax.add_feature(cfeature.LAKES)
  ax.add_feature(cfeature.COASTLINE)
  cax = ax.scatter(station_lon+360.0-central_longitude,station_lat,c=cluster_error,vmax=vmax,vmin=vmin,transform=crs)
  cb = fig.colorbar(cax)
  plt.title(method+' '+metric+' values for k='+str(k))
  fig.savefig('cluster_'+metric+'_'+method+'_k='+str(k)+'.png',bbox_inches='tight')
  plt.close()

  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.hist(cluster_error,range=(vmin,vmax),bins=20,cumulative=True)
  ax2 = ax.twinx()
  ax_ylim = ax.get_ylim()
  ax2.set_ylim([0,ax_ylim[1]/cluster_error.size*100])
  plt.title(method+' '+metric+' cdf for k='+str(k))
  fig.savefig('cdf_'+metric+'_'+method+'_k='+str(k)+'.png',bbox_inches='tight')
  plt.close()

  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.hist(cluster_error,range=(vmin,vmax),bins=20)
  ax2 = ax.twinx()
  ax_ylim = ax.get_ylim()
  ax2.set_ylim([0,ax_ylim[1]/cluster_error.size*100])
  plt.title(method+' '+metric+' histogram for k='+str(k))
  fig.savefig('hist_'+metric+'_'+method+'_k='+str(k)+'.png',bbox_inches='tight')
  plt.close()

##########################################################################
##########################################################################

def calculate_R2(x,y):
  
  xbar = np.mean(x)
  ybar = np.mean(y)
  top = np.dot(x-xbar,y-ybar)
  bottom = np.sqrt(np.dot(x-xbar,x-xbar))*np.sqrt(np.dot(y-ybar,y-ybar))
  r_value = top/bottom

  return r_value**2

##########################################################################
##########################################################################

def calculate_SMAPE(x,y):

  return np.mean(np.divide(np.absolute(x-y),np.absolute(x)+np.absolute(y)))

##########################################################################
##########################################################################
# Main program
##########################################################################
##########################################################################

if __name__ == '__main__':

  # Set name of pre-computed average file 
  if normalize:
    filename = 'averaged_month_'+str(month).zfill(2)+'_normalized_spectrum.nc'
  else:
    filename = 'averaged_month_'+str(month).zfill(2)+'_spectrum.nc'

  # Read in spectral data (perform averaging if pre-computed file doesn't exist)
  if os.path.isfile(filename):
    run_data,run_stations,model_spectra = read_averaged_model_spectra(filename)
  else:
    run_data,run_stations,model_spectra = get_averaged_model_spectra(run_direcs,month,normalize)
    write_averaged_model_spectra(filename,run_data,run_stations,model_spectra)
  
  for method in cluster_methods:
    for k in k_vals: 

      # Preform clustering
      labels,clustered_spectra,k = cluster_spectra(k,method,model_spectra,run_stations)
      
      # Compute colormap
      cmap = get_colormap(k)
    
      # Plot grid of cluster centers
      plot_cluster_grid(k,method,clustered_spectra,run_data['theta'],run_data['freq'],cmap)
       
      # Compute and plot the error between station spectra and cluster centers
      for metric in error_metrics:
        plot_cluster_error(k,method,metric,labels,model_spectra,clustered_spectra,run_stations['lon'],run_stations['lat'])  
    
      # Plot scatter plot of clustered stations
      plot_station_scatter(k,method,labels,run_stations['lon'],run_stations['lat'],run_data['theta'],run_data['freq'],clustered_spectra,model_spectra,cmap)



