import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.basemap import Basemap
from scipy import interpolate
import pprint
import glob
import netCDF4
import datetime
import os
plt.switch_backend('agg')

mode = 'validation'
#mode = 'stats'
#mode = 'individual'
year = '2002'
data_direc = '../../../../observed_data/spectrum/'
run_direc = '../../model_data/spectrum/'

variables = ['swden','swdir','swdir2','swr1','swr2']

#######################################################################
#######################################################################

def read_station_data(sta,year,variables):

  obs_data = {}
  success = True
  for i,var in enumerate(variables):

    # Open station data file
    filename = data_direc+sta+'_'+year+'_'+var+'.txt' 
    if not os.path.exists(filename):
      success = False
      print 'file not found: '+filename
      break 

    print '  reading '+var
    f = open(filename,'r')
    obs = f.read().splitlines()
  
    # Get the fequency data and initialize time list
    obs_data[var] = {}
    obs_data[var]['data'] = []
    if obs[0].find('mm') > 0:   # Adjust for different date formats
      n = 5                     # YYYY MM DD hh mm
    else:                       #      vs.
      n = 4                     # YYYY MM DD hh
    obs_data[var]['freq'] = obs[0].strip().split()[n:]
    obs_data[var]['datetime'] = []

    # Get the spectral variable data
    for j,line in enumerate(obs[1:]):
      line_sp = line.split()

      obs_data[var]['datetime'].append(' '.join(line_sp[0:n]))
      obs_data[var]['data'].append(line_sp[n:])
 
    # Convert observation data to numpy array
    try:
      obs_data[var]['data'] = np.asarray(obs_data[var]['data'],dtype=np.float)            
    except:
      print '  data shape not correct'
      success = False
    obs_data[var]['freq'] = np.asarray(obs_data[var]['freq'],dtype=np.float)
    
  
  return obs_data,success

#######################################################################
#######################################################################

def read_model_data(files,mode=None,stations=[]):
  
  for i,name in enumerate(files):
    print name.split('/')[-1]

    nc_file = netCDF4.Dataset(name,'r')

    # Initializations for first iteration of run
    if i == 0:
        # Station information
        if not stations:
          stations = {}
          stations['name'] = netCDF4.chartostring(nc_file.variables['station_name'][:,:]).tolist()
          stations['lon'] = np.squeeze(nc_file.variables['longitude'][:])
          stations['lat'] = np.squeeze(nc_file.variables['latitude'][:])
          nstations = len(stations['name'])
          station_ind = slice(nstations)
          stations2 = stations
        else:
          stations2 = {}
          model_stations = netCDF4.chartostring(nc_file.variables['station_name'][:,:]).tolist()
          station_ind = [model_stations.index(x) for x in stations['name'] if x in model_stations] 
          stations2['name'] = [x for x in stations['name'] if x in model_stations]
          nstations = len(station_ind)
          stations2['lon'] = np.squeeze(nc_file.variables['longitude'][0,station_ind])
          stations2['lat'] = np.squeeze(nc_file.variables['latitude'][0,station_ind])

        # Model output data (nfiles = number of timesnaps, 1 timesnap per file)
        nfiles = len(files)
        data = {}
        data['freq'] = nc_file.variables['frequency'][:] 
        data['theta'] = nc_file.variables['direction'][:]
        theta_ind = np.argsort(data['theta'])
        data['theta'] = np.radians(data['theta'][theta_ind])
        data['theta'] = np.insert(data['theta'],0,0.0)
        data['theta'] = np.append(data['theta'],2.0*np.pi)
        nfreq = len(data['freq'])
        ntheta = len(data['theta'])
        if mode == 'average':
          data['spectrum'] = np.zeros((1,nstations,nfreq,ntheta))
        else:
          data['spectrum'] = np.zeros((nfiles,nstations,nfreq,ntheta))

        # Time information
        data['time'] = np.empty(nfiles)
        data['ref_date'] = nc_file.variables['time'].getncattr('units').replace('days since ','')
        data['ref_date'] = datetime.datetime.strptime(data['ref_date'],'%Y-%m-%d %H:%M:%S')

    # Get time and output variables
    if i == 0:
      data['time'] = nc_file.variables['time'][:]
    efth = nc_file.variables['efth'][0,station_ind,:,theta_ind]
    nc_file.close()

    # Add averaged data for the  0,2*pi directions
    efth_0avg = 0.5*(efth[:,:,0] + efth[:,:,-1])
    efth = np.dstack((efth_0avg,efth))
    efth = np.dstack((efth,efth_0avg))
  
    if mode == 'average':
      data['spectrum'][0,:,:,:] = data['spectrum'][0,:,:,:] + efth 
    else:
      data['spectrum'][i,:,:,:] = efth 

  if mode == 'average':
    data['spectrum'] = data['spectrum']/float(nfiles)

  data['datetime'], data['date'] = output_time_to_date(data['time'],data['ref_date'])

  return data, stations2

################################################################################################
################################################################################################

def output_time_to_date(output_time,ref_date):

  output_date = []
  output_datetime = []

  # Convert output times to date format
  for i,t in enumerate(output_time):
    date = ref_date + datetime.timedelta(days=t) 

    output_date.append(date.strftime('%Y %m %d %H %M'))
    output_datetime.append(date)
  output_datetime = np.asarray(output_datetime,dtype='O')

  return output_datetime, output_date 

#######################################################################
#######################################################################

def interpolate_station_data(obs_data,variables):

  # Fix inconsistancies between number of frequency bins and data
  # (sometimes there are frequency bins listed with no data associated with them)
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
 
  # ----------------------------------------------------

  # Find coarsest frequency range to interpolate data to
  var_min = ''
  nfreq_min = 999 
  for var in variables:
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

def match_dates(obs_data,variables):

  success = True

  # Find the intersection of all dates across variables
  dates = set(obs_data['swden']['datetime'])
  for var in variables:
    dates = set.intersection(dates,set(obs_data[var]['datetime']))
  dates = list(dates)
  
  for var in variables:
    print '  matching dates '+var

    # Find indices of dates in intersection set
    indices = [i for i, x in enumerate(obs_data[var]['datetime']) if x in dates]
    indices = np.asarray(indices,dtype=np.int)
    
    # Select only dates that are in the intersection set 
    obs_data[var]['data'] = obs_data[var]['data'][indices,:]
    obs_data[var]['datetime'] = dates

  # Indicate whether any dates were matching 
  if len(dates) == 0:
    success = False

  return success

  
#######################################################################
#######################################################################

def compute_spectrum(obs_data,variables,mode='average'):

  print '  computing spectrum'

  # Setup spectrum dimensions
  nsnaps = len(obs_data['swden']['datetime'])
  nfreq = obs_data['swden']['freq'].size
  ndir = 37
  theta = np.radians(np.linspace(0.0,360.0,ndir))[np.newaxis].T
  obs_data['theta'] = theta

  # Initialize the spectrum variable
  if mode == 'average':
    spectrum = np.zeros((1,ndir,nfreq))
  elif mode == 'individual':
    spectrum = np.zeros((nsnaps,ndir,nfreq))
  elif mode == 'max':
    spectrum = np.zeros((1,ndir,nfreq))+-999.0
 
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
    elif mode == 'max':
      if not missing_data:
        spectrum[0,:,:] = np.maximum(spectrum[0,:,:],spec) 
      
    
  # Complete the average calculation
  if mode == 'average':
    spectrum = spectrum/float(nsnaps)

  return spectrum


#######################################################################
#######################################################################

def plot_station_spectrum(sta,lon,lat,lon2,lat2,freq1,theta1,spectrum1,labels1,freq2=[],theta2=[],spectrum2=None,labels2=None):

  print '  plotting spectrum'

  nsnaps = len(labels1)

  R1,Theta1 = np.meshgrid(freq1,theta1)
  if len(freq2) > 0 :
    R2,Theta2 = np.meshgrid(freq2,theta2)
  

  for i in range(nsnaps):

    if i % 20 == 0:

      # Create figure
      fig = plt.figure(figsize=[9.5,8])
      gs = gridspec.GridSpec(nrows=2,ncols=2,figure=fig)

      # Plot global station location
      ax = fig.add_subplot(gs[0,0])
      m = Basemap(projection='cyl',llcrnrlat= -90,urcrnrlat=90,\
                                   llcrnrlon=-180,urcrnrlon=180,resolution='c')
      m.fillcontinents(color='tan',lake_color='lightblue')
      m.drawcoastlines()
      ax.plot(lon,lat,'ro')
      ax.plot(lon2,lat2,'bo')

      # Plot local station location
      ax = fig.add_subplot(gs[0,1])
      m = Basemap(projection='cyl',llcrnrlat=lat-7.0 ,urcrnrlat=lat+7.0,\
                                   llcrnrlon=lon-10.0,urcrnrlon=lon+10.0,resolution='l')
      m.fillcontinents(color='tan',lake_color='lightblue')
      m.drawcoastlines()
      ax.plot(lon,lat,'ro')
      ax.plot(lon2,lat2,'bo')

      # Plot spectrum
      if spectrum2 is None:
        ax = fig.add_subplot(gs[1,:],polar=True)
      else:
        ax = fig.add_subplot(gs[1,0],polar=True)

      ax.set_theta_zero_location("N")
      ax.set_theta_direction(-1)
      cax = ax.contourf(Theta1,R1,spectrum1[i,:,:],30)
      cb = fig.colorbar(cax)
      ax.set_title('Directional spectrum '+labels1[i])
    
      if spectrum2 is not None :
     
        # Interploate spectrum onto data mesh with evenly spaced freq axis
        interp = interpolate.RegularGridInterpolator((theta2, freq2), spectrum2[i,:,:], bounds_error=False, fill_value=np.nan)
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
        ax = fig.add_subplot(gs[1,1],polar=True)
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        cax = ax.contourf(Theta1,R1,spec2,30)
        cb = fig.colorbar(cax)
        ax.set_title('Directional spectrum '+labels2[i])

      st = plt.suptitle('Station '+sta,fontsize=16)
      fig.tight_layout()
      fig.savefig(sta+'_'+year+'_'+str(i)+'.png',bbox_inches='tight',bbox_extra_artists=(st,))

      plt.close()

################################################################################################
################################################################################################

def read_station_file(station_file):

  # Read in stations names and location
  f = open(station_file)
  lines = f.read().splitlines()
  stations = {}
  stations['name'] = []
  stations['lon'] = []
  stations['lat'] = []
  for sta in lines:
    val = sta.split()
    stations['name'].append(val[2].strip("'"))
    stations['lon'].append(float(val[0]))
    stations['lat'].append(float(val[1]))
  nstations = len(stations['name'])
  stations['lon'] = np.asarray(stations['lon'])
  stations['lat'] = np.asarray(stations['lat'])

  return stations

#######################################################################
#######################################################################


if __name__ == '__main__':

  stations = read_station_file(data_direc+'stations.txt')

  run_files = sorted(glob.glob(run_direc+'spec*.nc'))
  if len(run_files) > 0:
    run_data,run_stations = read_model_data(run_files,'average',stations)
    model_spectrum_avg = np.zeros((1,run_data['theta'].size,run_data['freq'].size))


  for i,sta in enumerate(stations['name']):
  
    if sta.find(',') > 0:
      continue
  
    # Find station ID
    try:
      j = run_stations['name'].index(sta)
    except:
      continue

    print sta
    lon = stations['lon'][i]
    lat = stations['lat'][i]
    lon2 = run_stations['lon'][j]
    lat2 = run_stations['lat'][j]

    # Read station data variables
    obs_data,success = read_station_data(sta,year,variables)
    if not success:
      continue

    # Resolve differences in frequency bins across variables
    interpolate_station_data(obs_data,variables)
    
    # Resolve differences in obervation dates across variables 
    success = match_dates(obs_data,variables)

    if success:

      # Compute the spectrum from the data variables
      if mode == 'stats':
        spectrum_avg = compute_spectrum(obs_data,variables,'average')
        spectrum_max = compute_spectrum(obs_data,variables,'max')
      if mode == 'individual':
        spectrum = compute_spectrum(obs_data,variables,'individual')
      if mode == 'validation':
        obs_spectrum_avg = compute_spectrum(obs_data,variables,'average')
        


      # Plot the spectrum 
      if mode == 'stats':
        plot_station_spectrum(sta,lon,lat,obs_data['freq'],obs_data['theta'],spectrum_avg,['averaged over '+year], 
                                          obs_data['freq'],obs_data['theta'],spectrum_max,['maximum over '+year])  
      elif mode == 'individual':
        plot_station_spectrum(sta,lon,lat,obs_data['freq'],obs_data['theta'],spectrum,obs_data['datetime'])  
      elif mode == 'validation':
        model_spectrum_avg[0,:,:] = run_data['spectrum'][0,j,:,:].T
        plot_station_spectrum(sta,lon,lat,lon2,lat2,obs_data['freq'],obs_data['theta'],obs_spectrum_avg,['observations averaged over '+year],
                                          run_data['freq'],run_data['theta'],model_spectrum_avg,['model averaged over '+year])  
          
