import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
import matplotlib.ticker as ticker
from matplotlib import rc
import netCDF4 
import pprint
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import plot_alphabeta
import yaml
import os
plt.switch_backend('agg')



########################################################################
########################################################################

def write_alpha_beta(filename,header,nodes,lon,lat,sizes,alpha_spec,beta_spec):

  n = alpha_spec.shape[0]
  nfreq = alpha_spec.shape[1]
  ndir = alpha_spec.shape[2]

  lines = []
  lines.append(header)
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
  nfreq = alpha_spec.shape[1]
  ndir = alpha_spec.shape[2]

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

def plot_alpha_beta_spectra(Theta,Freq,alpha_spec,beta_spec,alpha_interp,beta_interp,kind):

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

  pwd = os.getcwd()

  f = open(pwd+'/rotate_alphabeta.config')
  cfg = yaml.load(f,yaml.Loader)
  pprint.pprint(cfg)
  
  theta = np.radians(np.linspace(0.0,360.0,cfg['ndir'],endpoint=False))
  freq = np.linspace(0.0,1.0,cfg['nfreq'])
  Theta,Freq = np.meshgrid(theta,freq)

  data = np.loadtxt(cfg['angled_file'])
  angled = data[:,2]

  lon_local,lat_local,nodes,sizes,alpha_local_avg,beta_local_avg,alpha_spec,beta_spec = plot_alphabeta.read_alpha_beta(cfg['filename_local_in'])
  alpha_interp,beta_interp = rotate_and_interpolate(Theta,nodes,angled,alpha_spec,beta_spec)
  header = '$WAVEWATCH III LOCAL OBSTRUCTIONS'
  write_alpha_beta(cfg['filename_local_out'],header,nodes,lon_local,lat_local,sizes,alpha_interp,beta_interp)
  plot_alpha_beta_spectra(Theta,Freq,alpha_spec,beta_spec,alpha_interp,beta_interp,'local')
  
  
  lon_shadow,lat_shadow,nodes,sizes,alpha_shadow_avg,beta_shadow_avg,alpha_spec,beta_spec = plot_alphabeta.read_alpha_beta(cfg['filename_shadow_in'])
  alpha_interp,beta_interp = rotate_and_interpolate(Theta,nodes,angled,alpha_spec,beta_spec)
  header = '$WAVEWATCH III SHADOW OBSTRUCTIONS'
  write_alpha_beta(cfg['filename_shadow_out'],header,nodes,lon_shadow,lat_shadow,sizes,alpha_interp,beta_interp)
  plot_alpha_beta_spectra(Theta,Freq,alpha_spec,beta_spec,alpha_interp,beta_interp,'shadow')
  
