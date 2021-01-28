import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
import matplotlib.ticker as ticker
from matplotlib import rc
import netCDF4 
import pprint
import cartopy.crs as ccrs
import cartopy.feature as cfeature
plt.switch_backend('agg')

rc('font',**{'family':'serif','serif':['Times New Roman']})

nfreq = 50

def read_alpha_beta(filename):
  f = open(filename,'r')
  lines = f.read().splitlines()
  
  n = int(lines[1])
  
  lon = []
  lat = []
  alpha = []
  beta = []
  line = 1
  for i in range(n):
   
    line = line + 1
   
    text = lines[line]
    text_sp = text.split()
    x = float(text_sp[7].replace(',',''))
    y = float(text_sp[9])
  
    line = line + 1 # node number
    line = line + 1 # sizes comment
    line = line + 1 # sizes
    line = line + 1 # mean alpha
  
    text = lines[line]
    text_sp = text.split()
    a = float(text_sp[-1])
  
    line = line + 1 # mean beta
    text = lines[line]
    text_sp = text.split()
    b = float(text_sp[-1])
  
    line = line + 1 # alpha comment
    line = line + 51 # beta comment
    line = line + 50 # lon lat comment
   
    lon.append(x)
    lat.append(y)
    alpha.append(a)
    beta.append(b)
  
  lon = np.array(lon)
  lat = np.array(lat)
  alpha = np.array(alpha)
  beta = np.array(beta)

  return [lon,lat,alpha,beta]

filename_local = 'obstructions_local.glo_unst.in'
filename_shadow = 'obstructions_shadow.glo_unst.in'
cmap = 'YlOrRd'
opac = 1.0 
land_color = '#DEB887'

lon,lat,alpha,beta = read_alpha_beta(filename_local)

fig = plt.figure(figsize=(6,4.5))
ax1 = fig.add_subplot(2,2,1,projection=ccrs.PlateCarree())
sc = ax1.scatter(lon,lat,c=alpha,vmin=0.0,vmax=1.0,cmap=cmap,transform=ccrs.PlateCarree(),alpha=opac,s=1)
ax1.add_feature(cfeature.OCEAN,zorder=0)
ax1.add_feature(cfeature.LAND, zorder=0, facecolor=land_color)
ax1.add_feature(cfeature.COASTLINE,alpha=0.5,zorder=0)
ax1.set_title('local alpha')

ax2 = fig.add_subplot(2,2,2,projection=ccrs.PlateCarree())
sc = ax2.scatter(lon,lat,c=beta,vmin=0.0,vmax=1.0,cmap=cmap,transform=ccrs.PlateCarree(),alpha=opac,s=1)
ax2.add_feature(cfeature.OCEAN,zorder=0)
ax2.add_feature(cfeature.LAND, zorder=0, facecolor=land_color)
ax2.add_feature(cfeature.COASTLINE,alpha=0.5,zorder=0)
ax2.set_title('local beta')


lon,lat,alpha,beta = read_alpha_beta(filename_shadow)

ax3 = fig.add_subplot(2,2,3,projection=ccrs.PlateCarree())
sc = ax3.scatter(lon,lat,c=alpha,vmin=0.0,vmax=1.0,cmap=cmap,transform=ccrs.PlateCarree(),alpha=opac,s=1)
ax3.add_feature(cfeature.OCEAN,zorder=0)
ax3.add_feature(cfeature.LAND, zorder=0, facecolor=land_color)
ax3.add_feature(cfeature.COASTLINE,alpha=0.5, zorder=0)
ax3.set_title('shadow alpha')

ax4 = fig.add_subplot(2,2,4,projection=ccrs.PlateCarree())
sc = ax4.scatter(lon,lat,c=beta,vmin=0.0,vmax=1.0,cmap=cmap,transform=ccrs.PlateCarree(),alpha=opac,s=1)
ax4.add_feature(cfeature.OCEAN,zorder=0)
ax4.add_feature(cfeature.LAND, zorder=0, facecolor=land_color)
ax4.add_feature(cfeature.COASTLINE,alpha=0.5, zorder=0)
ax4.set_title('shadow beta')

cb = fig.colorbar(sc,  ax=[ax1,ax2,ax3,ax4], orientation='horizontal',pad=0.05)
cb.set_label('alpha/beta value')

plt.savefig('alphabeta.png',bbox_inches='tight',dpi=400)
plt.close()
