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
ndir = 36
filename_local = 'obstructions_local.glo_unst.in'
filename_shadow = 'obstructions_shadow.glo_unst.in'
cmap = 'YlOrRd'
opac = 1.0 
land_color = '#DEB887'

theta = np.radians(np.linspace(0.0,360.0,ndir,endpoint=False))
freq = np.linspace(0.0,1.0,nfreq)
Theta,Freq = np.meshgrid(theta,freq)

#data = np.loadtxt('angled.d')
#angled = data[:,2]

########################################################################
########################################################################

def read_alpha_beta(filename):
  f = open(filename,'r')
  lines = f.read().splitlines()
  
 
  nodes = [] 
  lon = []
  lat = []
  sizes = []
  alpha_avg = []
  beta_avg = []
  alpha_spec = []
  beta_spec = []

  line = 1 # header comment
  n = int(lines[line])
  for i in range(n):
   
    line = line + 1 # lon lat comment
   
    text = lines[line]
    text_sp = text.split()
    x = float(text_sp[7].replace(',',''))
    y = float(text_sp[9])
  
    line = line + 1 # node number
    nodes.append(int(lines[line].split()[0]))

    line = line + 1 # sizes comment
    line = line + 1 # sizes
    sizes.append(lines[line])

    line = line + 1 # mean alpha
    text = lines[line]
    text_sp = text.split()
    a = float(text_sp[-1])
  
    line = line + 1 # mean beta
    text = lines[line]
    text_sp = text.split()
    b = float(text_sp[-1])
  
    line = line + 1 # alpha comment
    spectrum = []
    for i in range(nfreq):
      line = line + 1
      spectrum.append(lines[line].split())
    alpha_spec.append(spectrum)
    del(spectrum)

    line = line + 1 # beta comment
    spectrum = []
    for i in range(nfreq):
      line = line + 1
      spectrum.append(lines[line].split())
    beta_spec.append(spectrum)
    del(spectrum)
  
    lon.append(x)
    lat.append(y)
    alpha_avg.append(a)
    beta_avg.append(b)
  
  lon = np.array(lon)
  lat = np.array(lat)
  nodes = np.array(nodes)
  alpha_avg = np.array(alpha_avg)
  beta_avg = np.array(beta_avg)
  alpha_spec = np.array(alpha_spec)
  beta_spec = np.array(beta_spec)

  return [lon,lat,nodes,sizes,alpha_avg,beta_avg,alpha_spec,beta_spec]

########################################################################
########################################################################

def write_alpha_beta(filename,header,nodes,lon,lat,sizes,alpha_spec,beta_spec):

  lines = []
  lines.append(header)

  n = alpha_spec.shape[0]
  lines.append(str(n))

  for i in range(n):
    lines.append('$ ilon ilat of the cell. lon: {:.8f}, lat: {:.8f}'.format(lon[i],lat[i]))
    lines.append(str(nodes[i])+'   1')
    lines.append(sizes[i])
    lines.append('$ mean alpha: {:.16}'.format(np.mean(alpha_spec[i,:,:])))
    lines.append('$ mean beta: {:.16}'.format(np.mean(beta_spec[i,:,:])))
    lines.append('$alpha by ik, ith')
    for j in range(nfreq):
      line = ''
      for k in range(ndir):
        line = line + '{:.2f}  '.format(alpha_spec[i,j,k])
      lines.append(line)
    lines.append('$beta by ik, ith')
    for j in range(nfreq):
      line = ''
      for k in range(ndir):
        line = line + '{:.2f}  '.format(beta_spec[i,j,k])
      lines.append(line)
    
  f = open(filename,'w')
  for line in lines:
    f.write(line+'\n')
  f.close() 

########################################################################
########################################################################

def rotate_and_interpolate(Theta,nodes,angled,alpha_spec,beta_spec):

  n = alpha_spec.shape[0]
  alpha_interp = np.zeros((n,nfreq,ndir))
  beta_interp = np.zeros((n,nfreq,ndir))
  for i in range(n):
    nd = nodes[i]-1
    Theta2 = Theta - angled[nd]*np.pi/180.0
    for j in range(nfreq):
      alpha_interp[i,j,:] = np.interp(Theta2[j,:],Theta[j,:],alpha_spec[i,j,:],period = 2.0*np.pi)
      beta_interp[i,j,:]  = np.interp(Theta2[j,:],Theta[j,:],beta_spec[i,j,:], period = 2.0*np.pi)

  return [alpha_interp,beta_interp]

########################################################################
########################################################################

def plot_alpha_beta_spectra(alpha_spec,beta_spec,alpha_interp,beta_interp,kind):

  for i in range(10):
    print(i)
  
    fig = plt.figure(figsize=[8,4])
    ax = fig.add_subplot(2,2,1,polar=True)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.contourf(Theta,Freq,alpha_spec[i,:,:],30)
  
    ax = fig.add_subplot(2,2,2,polar=True)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.contourf(Theta,Freq,beta_spec[i,:,:],30)
  
    ax = fig.add_subplot(2,2,3,polar=True)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.contourf(Theta,Freq,alpha_interp[i,:,:],30)
    ax.set_title('AnglD = '+str(angled[nodes[i]-1]))
  
    ax = fig.add_subplot(2,2,4,polar=True)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.contourf(Theta,Freq,beta_interp[i,:,:],30)
    ax.set_title('AnglD = '+str(angled[nodes[i]-1]))
  
    plt.savefig(kind+'_spec_'+str(i)+'.png',bbox_inches='tight')
  
########################################################################
########################################################################

if __name__ == "__main__":


  lon_local,lat_local,nodes,sizes,alpha_local_avg,beta_local_avg,alpha_spec,beta_spec = read_alpha_beta(filename_local)
  #alpha_interp,beta_interp = rotate_and_interpolate(Theta,nodes,angled,alpha_spec,beta_spec)
  #filename = 'obstructions_local.glo_unst_RTDplus.in'
  #header = '$WAVEWATCH III LOCAL OBSTRUCTIONS'
  #write_alpha_beta(filename,header,nodes,lon_local,lat_local,sizes,alpha_interp,beta_interp)
  #plot_alpha_beta_spectra(alpha_spec,beta_spec,alpha_interp,beta_interp,'local')
  #
  #
  lon_shadow,lat_shadow,nodes,sizes,alpha_shadow_avg,beta_shadow_avg,alpha_spec,beta_spec = read_alpha_beta(filename_shadow)
  #alpha_interp,beta_interp = rotate_and_interpolate(Theta,nodes,angled,alpha_spec,beta_spec)
  #filename = 'obstructions_shadow.glo_unst_RTDplus.in'
  #header = '$WAVEWATCH III SHADOW OBSTRUCTIONS'
  #write_alpha_beta(filename,header,nodes,lon_shadow,lat_shadow,sizes,alpha_interp,beta_interp)
  #plot_alpha_beta_spectra(alpha_spec,beta_spec,alpha_interp,beta_interp,'shadow')
  
  
  fig = plt.figure(figsize=(6,4.5))
  ax1 = fig.add_subplot(2,2,1,projection=ccrs.PlateCarree())
  sc = ax1.scatter(lon_local,lat_local,c=alpha_local_avg,vmin=0.0,vmax=1.0,cmap=cmap,transform=ccrs.PlateCarree(),alpha=opac,s=1)
  ax1.add_feature(cfeature.OCEAN,zorder=0)
  ax1.add_feature(cfeature.LAND, zorder=0, facecolor=land_color)
  ax1.add_feature(cfeature.COASTLINE,alpha=0.5,zorder=0)
  ax1.set_title('local alpha')
  
  ax2 = fig.add_subplot(2,2,2,projection=ccrs.PlateCarree())
  sc = ax2.scatter(lon_local,lat_local,c=beta_local_avg,vmin=0.0,vmax=1.0,cmap=cmap,transform=ccrs.PlateCarree(),alpha=opac,s=1)
  ax2.add_feature(cfeature.OCEAN,zorder=0)
  ax2.add_feature(cfeature.LAND, zorder=0, facecolor=land_color)
  ax2.add_feature(cfeature.COASTLINE,alpha=0.5,zorder=0)
  ax2.set_title('local beta')
  
  ax3 = fig.add_subplot(2,2,3,projection=ccrs.PlateCarree())
  sc = ax3.scatter(lon_shadow,lat_shadow,c=alpha_shadow_avg,vmin=0.0,vmax=1.0,cmap=cmap,transform=ccrs.PlateCarree(),alpha=opac,s=1)
  ax3.add_feature(cfeature.OCEAN,zorder=0)
  ax3.add_feature(cfeature.LAND, zorder=0, facecolor=land_color)
  ax3.add_feature(cfeature.COASTLINE,alpha=0.5, zorder=0)
  ax3.set_title('shadow alpha')
  
  ax4 = fig.add_subplot(2,2,4,projection=ccrs.PlateCarree())
  sc = ax4.scatter(lon_shadow,lat_shadow,c=beta_shadow_avg,vmin=0.0,vmax=1.0,cmap=cmap,transform=ccrs.PlateCarree(),alpha=opac,s=1)
  ax4.add_feature(cfeature.OCEAN,zorder=0)
  ax4.add_feature(cfeature.LAND, zorder=0, facecolor=land_color)
  ax4.add_feature(cfeature.COASTLINE,alpha=0.5, zorder=0)
  ax4.set_title('shadow beta')
  
  cb = fig.colorbar(sc,  ax=[ax1,ax2,ax3,ax4], orientation='horizontal',pad=0.05)
  cb.set_label('alpha/beta value')
  
  plt.savefig('alphabeta.png',bbox_inches='tight',dpi=400)
  plt.close()
