import numpy as np
import netCDF4
from scipy import interpolate
import matplotlib.pyplot as plt
import datetime
import subprocess
import os
import yaml
import pprint

plt.switch_backend('agg')

########################################################################
# Read in data
########################################################################

pwd = os.getcwd()
f = open(pwd+'/prep_icebergs.config')
cfg = yaml.load(f)
pprint.pprint(cfg)

# Get data from ice file
ice_nc = netCDF4.Dataset(cfg['ice_file'],'r')
lon_ice = ice_nc.variables['lon'][:]
lat_ice = ice_nc.variables['lat'][:]
time_ice = ice_nc.variables['time'][:].astype(np.float64)
shape_ice = (time_ice.size, lat_ice.size, lon_ice.size)

if cfg['iceberg_file'] != '':
  # Get data from iceberg file
  iceberg_nc = netCDF4.Dataset(cfg['iceberg_file'],'r')
  lon_iceberg = iceberg_nc.variables['longitude'][:]
  lat_iceberg = iceberg_nc.variables['latitude'][:]
  time_iceberg = iceberg_nc.variables['time'][:].astype(np.float64)
  iceberg_probability = iceberg_nc.variables['probability'][:,:,:]
  iceberg_area = iceberg_nc.variables['ice_area'][:,:,:]
  
  # Extend iceberg data to first snap of the next year
  iceberg_nc_next = netCDF4.Dataset(cfg['iceberg_file_next'],'r')
  time_iceberg = np.append(time_iceberg,iceberg_nc_next.variables['time'][0].astype(np.float64))
  iceberg_probability = np.concatenate((iceberg_probability,np.expand_dims(iceberg_nc_next.variables['probability'][0,:,:],axis=0)))
  iceberg_area = np.concatenate((iceberg_area,np.expand_dims(iceberg_nc_next.variables['ice_area'][0,:,:],axis=0)))

  # Calculate area of grid cells from iceberg dataset  
  Rearth = 6371.229                                                                     
  deg2rad = np.pi/180.0                                                                 
  dLon = (lon_iceberg[1] - lon_iceberg[0])*deg2rad                                      
  lat_iceberg_cell = lat_iceberg + 0.5*(lat_iceberg[1]-lat_iceberg[0])                  
  lat_iceberg_cell = np.insert(lat_iceberg_cell,0,-90.1)                                
  dLat = np.sin(lat_iceberg_cell[1:]*deg2rad) - np.sin(lat_iceberg_cell[0:-1]*deg2rad)  
  Asw = Rearth**2*dLon*dLat                                                             

  # Calculate iceberg damping 
  W = 0.42                                                                           # effective iceberg width
  iceberg_damping = np.divide(np.multiply(iceberg_probability,iceberg_area),Asw)/W   # inverse of e-folding scale from Ardhuin et el. 2011

else:

  iceberg_damping = np.zeros(shape_ice)

########################################################################
# Spatial interpolation 
########################################################################

if iceberg_damping.shape != shape_ice:
  # Create high-res grid of points to interpolate to 
  Lon_ice,Lat_ice = np.meshgrid(lon_ice,lat_ice)
  shape = (time_iceberg.size, Lon_ice.shape[0], Lon_ice.shape[1])
  iceberg_damping_interp = np.zeros(shape)
  
  # Find sub-section of high-res grid that corresponds with the iceberg data
  min_lon_iceberg = np.amin(lon_iceberg)
  max_lon_iceberg = np.amax(lon_iceberg)
  min_lat_iceberg = np.amin(lat_iceberg)
  max_lat_iceberg = np.amax(lat_iceberg)
  idx_lon, = np.where((lon_ice >= min_lon_iceberg) & (lon_ice <= max_lon_iceberg))
  idx_lat, = np.where((lat_ice >= min_lat_iceberg) & (lat_ice <= max_lat_iceberg))
  subsection_idx = np.ix_(idx_lat,idx_lon)
  subsection_shape = Lon_ice[subsection_idx].shape
  
  # Get points to interpolate 
  Lon = Lon_ice[subsection_idx].ravel()
  Lat = Lat_ice[subsection_idx].ravel()
  pts = np.vstack((Lon,Lat)).T
  
  # Interpolate iceberg timesnaps in space
  for i,t in enumerate(time_iceberg):
   
    # Replace nan values with average of non-nan neighbor values
    iceberg_damping_filled = np.copy(iceberg_damping[i,:,:])
    idx,idy = np.where(np.isnan(iceberg_damping[i,:,:]))
    for j in range(len(idx)):
      neighborhood = iceberg_damping[i,idx[j]-1:idx[j]+2,idy[j]-1:idy[j]+2]
      n = np.count_nonzero(~np.isnan(neighborhood))
      iceberg_damping_filled[idx[j],idy[j]] = np.nansum(neighborhood)/n
  
    # Interpolate and insert values back into high-res grid
    interpolate_icebergs = interpolate.RegularGridInterpolator((lon_iceberg,lat_iceberg),iceberg_damping_filled)
    interp_data = interpolate_icebergs(pts)
    iceberg_damping_interp[i][subsection_idx] = np.reshape(interp_data,subsection_shape)
    
    # Plot coarse grid iceberg concentration
    plt.figure(figsize=(12,4))
    plt.contourf(lon_iceberg,lat_iceberg,iceberg_damping[i,:,:].T)
    plt.ylim(-90.0,-40.0)
    plt.savefig('raw_icebergs_'+str(i).zfill(3)+'.png')
    plt.close()
    
    # Plot interpolated iceberg concentration
    plt.figure(figsize=(12,4))
    plt.contourf(lon_ice,lat_ice,iceberg_damping_interp[i,:,:])
    plt.ylim(-90.0,-40.0)
    plt.savefig('interp_icebergs_'+str(i).zfill(3)+'.png')
    plt.close()

else:
  iceberg_damping_interp = iceberg_damping

########################################################################
# Temporal interpolation 
########################################################################

if iceberg_damping_interp.shape != shape_ice:
  # Get reference times for ice and iceberg data
  ref_date_ice = ice_nc.variables['time'].getncattr('units').replace('hours since ','')
  ref_date_ice = datetime.datetime.strptime(ref_date_ice.replace('.0 +0:00',''),'%Y-%m-%d %H:%M:%S')
  ref_date_iceberg = iceberg_nc.variables['time'].getncattr('units').replace('days since ','')
  ref_date_iceberg = datetime.datetime.strptime(ref_date_iceberg,'%Y-%m-%dT%H:%M:%SZ')
  
  # Initialize iceberg time window
  i_iceberg = 0
  t_iceberg1 = ref_date_iceberg + datetime.timedelta(days=time_iceberg[i_iceberg])
  t_iceberg2 = ref_date_iceberg + datetime.timedelta(days=time_iceberg[i_iceberg+1])
  
  # Initialize final iceberg data array
  iceberg_damping_final = np.zeros(shape_ice)
  
  for i,t in enumerate(time_ice):
  
    # Find iceberg time window that corresponds with ice time
    t_ice = ref_date_ice + datetime.timedelta(hours=t)
    print datetime.datetime.strftime(t_ice,'%Y-%m-%d %H:%M:%S') 
    print datetime.datetime.strftime(t_iceberg1,'%Y-%m-%d %H:%M:%S'),datetime.datetime.strftime(t_iceberg2,'%Y-%m-%d %H:%M:%S')
    if (t_ice > t_iceberg2):
      i_iceberg = i_iceberg + 1
      t_iceberg1 = t_iceberg2
      t_iceberg2 = ref_date_iceberg + datetime.timedelta(days=time_iceberg[i_iceberg+1])
  
    # Interpolate in time
    if cfg['interp_type'] == 'linear':
      top = t_ice-t_iceberg2
      bottom = t_iceberg1-t_iceberg2
      phi1 = top.total_seconds()/bottom.total_seconds()
      top = t_ice-t_iceberg1
      bottom = t_iceberg2-t_iceberg1
      phi2 = top.total_seconds()/bottom.total_seconds()
      iceberg_damping_final[i,:,:] = phi1*iceberg_damping_interp[i_iceberg,:,:] + phi2*iceberg_damping_interp[i_iceberg+1,:,:]
    elif cfg['interp_type'] == 'constant':
      iceberg_damping_final[i,:,:] = iceberg_damping_interp[i_iceberg,:,:]
  
    # Plot timesnaps
    if i % 20 == 0:
      plt.figure()
      plt.contourf(lon_ice,lat_ice,ice_nc.variables['ICEC_L1'][i,:,:])
      plt.savefig('ice_'+str(i).zfill(4)+'.png')
      plt.close()
  
      plt.figure(figsize=(12,4))
      plt.contourf(lon_ice,lat_ice,iceberg_damping_final[i,:,:])
      plt.ylim(-90.0,-40.0)
      plt.savefig('final_icebergs_'+str(i).zfill(4)+'.png')
      plt.close()

else:

  iceberg_damping_final = iceberg_damping_interp
 
########################################################################
# Create combined NetCDF file 
########################################################################

# Copy ice file and open for appending
dest = pwd+'/'+cfg['combined_file']
subprocess.call(['cp',cfg['ice_file'],dest])
ncfile = netCDF4.Dataset(dest,'a')

# Add the iceberg data
if 'iceberg_damping' not in ncfile.variables.keys():
  nc_iceberg_damping = ncfile.createVariable('iceberg_damping','f8',('time','lat','lon'))
nc_iceberg_damping[:] = iceberg_damping_final
nc_iceberg_damping.units = 'km-1'
nc_iceberg_damping.long_name = 'iceberg damping'
nc_iceberg_damping.product_description = 'Portion of the incoming wave energy blocked by icebergs over a unit propagation distance. Calculated using Eqn (2) in Ardhuin et al. Observation and parameterization of small icebergs: Drifting breakwaters in the southern ocean. Ocean Modelling. 2011.'
ncfile.close()
