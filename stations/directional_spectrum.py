import numpy as np
import matplotlib.pyplot as plt
import pprint

data_direc = 'spectral_data/'
year = '2005'

variables = ['swden','swdir','swdir2','swr1','swr2']

#######################################################################
#######################################################################

def read_station_data(sta,year,variables):

  obs_data = {}
  for i,var in enumerate(variables):
    print '  '+var

    # Open station data file
    f = open(data_direc+sta+'_'+year+'_'+var+'.txt','r')
    obs = f.read().splitlines()
  
    # Get the fequency data and initialize time list
    obs_data[var] = []
    if i == 0:
      obs_data['freq'] = obs[0].split()[5:]
      obs_data['freq'] = np.asarray(obs_data['freq'])
      obs_data['freq'] = obs_data['freq'].astype(np.float)
      obs_data['datetime'] = []
      nfreq = len(obs_data['freq'])

    # Get the spectral variable data
    for j,line in enumerate(obs[1:]):
      line_sp = line.split()

      if i == 0:
        obs_data['datetime'].append(' '.join(line_sp[0:5]))
      obs_data[var].append(line_sp[5:])
      if len(line_sp[5:]) != nfreq:
        print '  dimension error in line '+str(j)
  
    # Convert observation data
    obs_data[var] = np.asarray(obs_data[var],dtype=np.float)            
    
  
  return obs_data

#######################################################################
#######################################################################

def plot_station_spectrum(sta,obs_data):


  nsnaps = len(obs_data['datetime'])

  theta = np.radians(np.linspace(0.0,360.0,37))[np.newaxis].T
  print theta.shape, obs_data['freq'].shape
  n = obs_data['freq'].size
  R,Theta = np.meshgrid(obs_data['freq'],theta)

  for i in range(nsnaps):
    C11    = obs_data['swden'][i,0:n]
    alpha1 = np.radians(obs_data['swdir'][i,0:n])
    alpha2 = np.radians(obs_data['swdir2'][i,0:n])
    r1     = obs_data['swr1'][i,0:n]*0.01
    r2     = obs_data['swr2'][i,0:n]*0.01

    spectrum = C11*(0.5 + r1*np.cos(theta-alpha1) + r2*np.cos(2.0*(theta-alpha2)))/np.pi

    fig,ax = plt.subplots(subplot_kw=dict(projection='polar'))
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    cax = ax.contourf(Theta,R,spectrum,30)
    cb = fig.colorbar(cax)
    fig.tight_layout()
    fig.savefig(sta+'_'+str(i)+'.png',bbox_inches='tight')

    if i > 15:
      break


#######################################################################
#######################################################################

if __name__ == '__main__':

  # Read the list of stations
  f = open(data_direc+'stations.txt','r')
  stations = f.read().splitlines()
  f.close()
  
  for i,line in enumerate(stations):
  

    # Find station ID
    sta = line.split('  ')[2].replace("'","")
    if sta != '42001':
      continue 

    print sta
  
    obs_data = read_station_data(sta,year,variables)
  
    plot_station_spectrum(sta,obs_data)  
      
    break
