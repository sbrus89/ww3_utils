import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.basemap import Basemap
import pprint

mode = 'average'
#mode = 'individual'
data_direc = 'spectral_data/'
year = '2005'

variables = ['swden','swdir','swdir2','swr1','swr2']

#######################################################################
#######################################################################

def read_station_data(sta,year,variables):

  obs_data = {}
  for i,var in enumerate(variables):
    print '  reading'+var

    success = True
 
    # Open station data file
    f = open(data_direc+sta+'_'+year+'_'+var+'.txt','r')
    obs = f.read().splitlines()
  
    # Get the fequency data and initialize time list
    obs_data[var] = []
    if i == 0:
      obs_data['freq'] = obs[0].strip().split()[5:]
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
        success = False
 
    # Convert observation data to numpy array
    obs_data[var] = np.asarray(obs_data[var],dtype=np.float)            
    
  
  return obs_data,success

#######################################################################
#######################################################################

def compute_spectrum(obs_data, mode='average'):

  print '  computing spectrum'

  success = True

  # Check that number of time snaps are the same for all data variables
  nsnaps = len(obs_data['datetime'])
  if obs_data['swden'].shape[0] != nsnaps or obs_data['swdir'].shape[0] != nsnaps or obs_data['swdir2'].shape[0] != nsnaps or obs_data['swr1'].shape[0] != nsnaps or obs_data['swr2'].shape[0] != nsnaps:
    success = False

  # Setup spectrum dimensions 
  nfreq = obs_data['freq'].size
  ndir = 37
  theta = np.radians(np.linspace(0.0,360.0,ndir))[np.newaxis].T
  obs_data['theta'] = theta

  # Initialize the spectrum variable
  if mode == 'average':
    spectrum = np.zeros((1,ndir,nfreq))
  elif mode == 'individual':
    spectrum = np.zeros((nsnaps,ndir,nfreq))
 
  # Early return if number of time snaps isn't the same for all data variables
  if success == False:
    return spectrum, success

  # Compute the spectrum for each time snap
  for i in range(nsnaps):
 
    # Get data variables
    C11    = obs_data['swden'][i,0:nfreq]
    alpha1 = np.radians(obs_data['swdir'][i,0:nfreq])
    alpha2 = np.radians(obs_data['swdir2'][i,0:nfreq])
    r1     = obs_data['swr1'][i,0:nfreq]*0.01
    r2     = obs_data['swr2'][i,0:nfreq]*0.01
    
    # Check for missing data
    missing_data = False
    if np.amax(C11) == 999.0 or np.amax(alpha1) == 999.0 or np.amax(alpha2) == 999.0 or np.amax(r1) == 999.0 or np.amax(r1) == 999.0:
      missing_data = True
 
    # Compute spectrum from data
    spec = C11*(0.5 + r1*np.cos(theta-alpha1) + r2*np.cos(2.0*(theta-alpha2)))/np.pi

    # Either store spectrum of accumulate for average
    if mode == 'average':
      if not missing_data:
        spectrum[0,:,:] = spectrum[0,:,:] + spec 
    elif mode == 'individual':
      spectrum[i,:,:] = spec 
    
  # Complete the average calculation
  if mode == 'average':
    spectrum = spectrum/float(nsnaps)

  return spectrum, success


#######################################################################
#######################################################################

def plot_station_spectrum(sta,lon,lat,dates,freq,theta,obs_data):

  print '  plotting spectrum'

  nsnaps = len(dates)

  R,Theta = np.meshgrid(freq,theta)

  for i in range(nsnaps):

    if i % 20 == 0:

      # Create figure
      fig = plt.figure(figsize=[6,8])
      gs = gridspec.GridSpec(nrows=2,ncols=2,figure=fig)

      # Plot global station location
      ax = fig.add_subplot(gs[0,0])
      m = Basemap(projection='cyl',llcrnrlat= -90,urcrnrlat=90,\
                                   llcrnrlon=-180,urcrnrlon=180,resolution='c')
      m.fillcontinents(color='tan',lake_color='lightblue')
      m.drawcoastlines()
      ax.plot(lon,lat,'ro')

      # Plot local station location
      ax = fig.add_subplot(gs[0,1])
      m = Basemap(projection='cyl',llcrnrlat=lat-7.0 ,urcrnrlat=lat+7.0,\
                                   llcrnrlon=lon-10.0,urcrnrlon=lon+10.0,resolution='l')
      m.fillcontinents(color='tan',lake_color='lightblue')
      m.drawcoastlines()
      ax.plot(lon,lat,'ro')

      # Plot spectrum
      ax = fig.add_subplot(gs[1,:],polar=True)
      ax.set_theta_zero_location("N")
      ax.set_theta_direction(-1)
      cax = ax.contourf(Theta,R,spectrum[i,:,:],30)
      cb = fig.colorbar(cax)
      ax.set_title('Directional spectrum '+dates[i])
      st = plt.suptitle('Station '+sta,fontsize=16)
      fig.tight_layout()
      fig.savefig(sta+'_'+str(i)+'.png',bbox_inches='tight',bbox_extra_artists=(st,))

      plt.close()


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
    lon = float(line.split('  ')[0])
    lat = float(line.split('  ')[1])
    print sta

    #if sta != '46042':
    #  continue
  
    # Read station data variables
    obs_data,success = read_station_data(sta,year,variables)
    
    # Compute the spectrum from the data variables
    if success:
      spectrum, success = compute_spectrum(obs_data,mode)
  
    # Plot the spectrum 
    if success:
      if mode == 'average':
        plot_station_spectrum(sta,lon,lat,['averaged over '+year],obs_data['freq'],obs_data['theta'],spectrum)  
      elif mode == 'individual':
        plot_station_spectrum(sta,lon,lat,obs_data['datetime'],obs_data['freq'],obs_data['theta'],spectrum)  
        
