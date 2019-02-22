import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.basemap import Basemap
from scipy import interpolate
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
    print '  reading '+var

    # Open station data file
    f = open(data_direc+sta+'_'+year+'_'+var+'.txt','r')
    obs = f.read().splitlines()
  
    # Get the fequency data and initialize time list
    obs_data[var] = {}
    obs_data[var]['data'] = []
    obs_data[var]['freq'] = obs[0].strip().split()[5:]
    obs_data[var]['datetime'] = []

    # Get the spectral variable data
    for j,line in enumerate(obs[1:]):
      line_sp = line.split()

      obs_data[var]['datetime'].append(' '.join(line_sp[0:5]))
      obs_data[var]['data'].append(line_sp[5:])
 
    # Convert observation data to numpy array
    obs_data[var]['data'] = np.asarray(obs_data[var]['data'],dtype=np.float)            
    obs_data[var]['freq'] = np.asarray(obs_data[var]['freq'],dtype=np.float)
    
  
  return obs_data

#######################################################################
#######################################################################

def interpolate_station_data(obs_data,variables):

  # Fix inconsistancies between number of frequency bins and data
  for var in variables:
    nfreq = obs_data[var]['data'].shape[1]
    obs_data[var]['freq'] = obs_data[var]['freq'][0:nfreq]

  # Check if number of frequency bins is different accross data variables
  nfreq = obs_data['swden']['freq'].size
  need_to_interpolate = False
  for var in variables:
    if obs_data[var]['freq'].size != nfreq:
      need_to_interpolate = True

  # Early return if interpolation is not needed
  if need_to_interpolate == False:
    obs_data['freq'] = obs_data['swden']['freq']
    return 
 

  # Find coarsest frequency range to interpolate data to
  var_min = ''
  nfreq_min = 999 
  for var in variables:
    print var, nfreq
    if nfreq < nfreq_min:
      var_min = var
      nfreq_min = nfreq

  # Find min and max of frequency ranges 
  freq_min = 999.0
  freq_max = -999.0
  for var in variables:
    if var != var_min:
      fmin = np.amin(obs_data[var]['freq'])
      fmax = np.amax(obs_data[var]['freq'])
      if fmin < freq_min:
        freq_min = fmin
      if fmax > freq_max:
        freq_max = fmax

  # Make sure frequency range used for interpolation is within other frequency ranges
  data_mask, = np.where((obs_data[var_min]['freq'] >= freq_min) &  (obs_data[var_min]['freq'] <= freq_max))
  obs_data['freq'] = obs_data[var_min]['freq'][data_mask]
  print obs_data['freq']

  # Interploate variables
  for var in variables:    
      print '  interpolating '+ var

      # Handle degree interpolation
      if var.find('dir') > 0:
        data = np.rad2deg(np.unwrap(np.deg2rad(obs_data[var]['data'])))
      else:
        data = obs_data[var]['data']
      
      # Create interpolation object
      f = interpolate.interp1d(obs_data[var]['freq'],data)     

      # Perform interpolation
      interp_data = f(obs_data['freq'])

      # Handle degree interpolation
      if var.find('dir') > 0:
        interp_data = interp_data % 360

      # Store interpolated data
      obs_data[var]['data'] = interp_data
      obs_data[var]['freq'] = obs_data['freq']

#######################################################################
#######################################################################

def compute_spectrum(obs_data,variables,mode='average'):

  print '  computing spectrum'

  success = True

  # Check that number of time snaps are the same for all data variables
  nsnaps = len(obs_data['swden']['datetime'])
  for var in variables:
    if obs_data[var]['data'].shape[0] != nsnaps:
      success = False

  # Setup spectrum dimensions 
  nfreq = obs_data['swden']['freq'].size
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
 
    # Check for missing data
    missing_data = False
    for var in variables:
      if np.amax(obs_data[var]['data'][i,0:nfreq]) == 999.0:
        missing_data = True

    # Get data variables
    C11    = obs_data['swden']['data'][i,0:nfreq]
    alpha1 = np.radians(obs_data['swdir']['data'][i,0:nfreq])
    alpha2 = np.radians(obs_data['swdir2']['data'][i,0:nfreq])
    r1     = obs_data['swr1']['data'][i,0:nfreq]*0.01
    r2     = obs_data['swr2']['data'][i,0:nfreq]*0.01
 
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
    obs_data = read_station_data(sta,year,variables)

    interpolate_station_data(obs_data,variables)
    
    # Compute the spectrum from the data variables
    spectrum, success = compute_spectrum(obs_data,variables,mode)
  
    # Plot the spectrum 
    if success:
      if mode == 'average':
        plot_station_spectrum(sta,lon,lat,['averaged over '+year],obs_data['freq'],obs_data['theta'],spectrum)  
      elif mode == 'individual':
        plot_station_spectrum(sta,lon,lat,obs_data['datetime'],obs_data['freq'],obs_data['theta'],spectrum)  
        
