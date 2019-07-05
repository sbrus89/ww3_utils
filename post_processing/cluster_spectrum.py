import glob
import sys
import numpy as np
from matplotlib import cm
from matplotlib.transforms import Bbox
from matplotlib.patches import Rectangle
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from scipy import interpolate
import directional_spectrum
from sklearn.cluster import KMeans
from sklearn.cluster import DBSCAN 
from sklearn.cluster import AgglomerativeClustering
from sklearn.mixture import GaussianMixture
from sklearn.decomposition import PCA
plt.switch_backend('agg')
np.set_printoptions(threshold=sys.maxsize)

stations_direc = '../../../../../'
run_direc = '../../model_data/spectrum/'

method = 'mixture'


if __name__ == '__main__':

  k = 25 


  stations = directional_spectrum.read_station_file(stations_direc+'stations_grid.txt')

  run_files = sorted(glob.glob(run_direc+'spec.200201*T00Z_spec.nc'))
  if len(run_files) > 0:
    #run_data,run_stations = directional_spectrum.read_model_data(run_files,'average',stations)
    run_data,run_stations = directional_spectrum.read_model_data(run_files,'average')

  print len(stations['name'])
  print len(run_stations['name'])
  one,nstations,nfrequency,ndirection =  run_data['spectrum'].shape 


  spec = np.reshape(run_data['spectrum'],(nstations,nfrequency*ndirection))
  print spec.shape
  
  if method == 'kmeans':
    clustered = KMeans(n_clusters=k,random_state=0).fit(spec)
    spec = clustered.cluster_centers_
    labels = clustered.labels_
  elif method == 'pca':
    clustered = PCA(n_components=k).fit(spec)
    spec = clustered.components_
  elif method == 'dbscan':
    clustered = DBSCAN(eps=1e-1,metric='l2',min_samples=100).fit(spec)
    spec = clustered.components_
    labels = clusterd.labels_
  elif method == 'mixture':
    clustered = GaussianMixture(n_components=k).fit(spec)
    spec = clustered.means_
    labels = clustered.predict(spec)

    
  print spec.shape
  k = spec.shape[0]

  cmap1 = cm.get_cmap('tab20b')(range(20))
  cmap2 = cm.get_cmap('tab20c')(range(20))
  cmap = ListedColormap(np.vstack((cmap1,cmap2)))
  cmap = ListedColormap(cmap(range(k)))

  if method == 'kmeans' or method == 'dbscan':

    fig = plt.figure(figsize=[18.0,6.0])
    ax = fig.add_subplot()
    m = Basemap(projection='cyl',llcrnrlat= -90,urcrnrlat=90,\
                                 llcrnrlon=-180,urcrnrlon=180,resolution='c')
    m.fillcontinents(color='tan',lake_color='lightblue')
    m.drawcoastlines()
    plt.scatter(run_stations['lon'],run_stations['lat'],c=labels,cmap=cmap)
    fig.tight_layout()
    fig.savefig('clustered_spectrum.png',bbox_inches='tight')
    plt.close()

    fig = plt.figure()
    x = np.arange(k)
    plt.scatter(x,x,c=x,cmap=cmap)
    fig.savefig('colors.png',bbox_inches='tight')
    plt.close()


  theta2 = run_data['theta']
  freq2 = run_data['freq']  

  ndir = 37
  nfreq = 48
  theta1 = np.radians(np.linspace(0.0,360.0,ndir))[np.newaxis].T
  freq1 = np.linspace(0.0,0.5,nfreq)
  R1,Theta1 = np.meshgrid(freq1,theta1)


  n = int(np.ceil(np.sqrt(k)))

  fig = plt.figure(figsize=[18,18])
  for i in range(k):
    spectrum = np.reshape(spec[i,:],(nfrequency,ndirection)).T  

    # Interploate spectrum onto data mesh with evenly spaced freq axis
    interp = interpolate.RegularGridInterpolator((theta2,freq2), spectrum, bounds_error=False, fill_value=np.nan)
    pts = np.vstack((Theta1.ravel(),R1.ravel())).T
    spec2 = interp(pts)
    spec2 = np.reshape(spec2,Theta1.shape)

    # Adjust for WW3 theta convention
    theta_flat = theta1[0:-1].flatten()
    cols =  np.asarray(range(0,len(freq1)))
    flipped = np.mod(-(theta_flat-np.pi/4.0)+np.pi/4.0,2.0*np.pi)         # Reflect across
    idx = np.argsort(flipped)                                             # the 45/225 degree line
    idx = np.ix_(idx,cols)                                                # and rearrange 
    spec2 = spec2[idx]
    flipped = np.mod(-(theta_flat-3.0*np.pi/4.0)+3.0*np.pi/4.0,2.0*np.pi) # Reflect across the
    idx = np.argsort(flipped)                                             # 135/315 degree line
    idx = np.ix_(idx,cols)                                                # and rearrange
    spec2 = spec2[idx]
    first_col =  spec2[0,:]
    spec2 = np.vstack((spec2,first_col[np.newaxis,:]))

    # Plot spectrum
    ax = fig.add_subplot(n,n,i+1,polar=True)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    cax = ax.contourf(Theta1,R1,spec2,30)
    cb = fig.colorbar(cax)


  fig.tight_layout()
  h = 1.0/float(n)
  c = 0
  for j in range(n):
    for i in range(n):
     rect = [Rectangle([0.0+(i*h),1.0-(j*h)-h],h,h,facecolor=cmap(c),edgecolor='none',zorder=-1,transform=fig.transFigure)]
     c = c +1
     fig.patches.extend(rect)
  
  fig.savefig('k_spectra.png',bbox_inches='tight')
