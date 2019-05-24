import numpy as np
import netCDF4
from scipy import interpolate
import matplotlib.pyplot as plt
plt.switch_backend('agg')

# Data file paths, ice is assumed to be higher resolution in space  
ice_file = '../ice/icecon.200201-200212.ww3.nc'
iceberg_file = './prod_latlon_jason1_2002.nc'

# Get data from ice file
ice_nc = netCDF4.Dataset(ice_file,'r')
lon_ice = ice_nc.variables['lon'][:]
lat_ice = ice_nc.variables['lat'][:]
time_ice = ice_nc.variables['time'][:]

# Get data from iceberg file
iceberg_nc = netCDF4.Dataset(iceberg_file,'r')
lon_iceberg = iceberg_nc.variables['longitude'][:]
lat_iceberg = iceberg_nc.variables['latitude'][:]
time_iceberg = iceberg_nc.variables['time'][:]

# Calculate iceberg concentration
iceberg_probability = iceberg_nc.variables['probability'][:,:,:]
iceberg_area = iceberg_nc.variables['ice_area'][:,:,:]
Asw = 100.0
W = 0.42
iceberg_concentration = np.multiply(iceberg_probability,iceberg_area)/Asw/W

# Create high-res grid of points to interpolate to 
Lon_ice,Lat_ice = np.meshgrid(lon_ice,lat_ice)
iceberg_concentration_interp = np.zeros(Lon_ice.shape)

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
  idx,idy = np.where(np.isnan(iceberg_concentration[i,:,:]))
  for j in range(len(idx)):
    neighborhood = iceberg_concentration[i,idx[j]-1:idx[j]+2,idy[j]-1:idy[j]+2]
    n = np.count_nonzero(~np.isnan(neighborhood))
    iceberg_concentration[i,idx[j],idy[j]] = np.nansum(neighborhood)/n

  # Interpolate and insert values back into high-res grid
  interpolate_icebergs = interpolate.RegularGridInterpolator((lon_iceberg,lat_iceberg),iceberg_concentration[i,:,:])
  interp_data = interpolate_icebergs(pts)
  iceberg_concentration_interp[subsection_idx] = np.reshape(interp_data,subsection_shape)
  
  # Plot coarse grid iceberg concentration
  plt.figure(figsize=(12,4))
  plt.contourf(lon_iceberg,lat_iceberg,iceberg_concentration[i,:,:].T)
  plt.ylim(-90.0,-40.0)
  plt.savefig('raw_icebergs_'+str(i).zfill(3)+'.png')
  plt.close()
  
  # Plot interpolated iceberg concentration
  plt.figure(figsize=(12,4))
  plt.contourf(lon_ice,lat_ice,iceberg_concentration_interp)
  plt.ylim(-90.0,-40.0)
  plt.savefig('interp_icebergs_'+str(i).zfill(3)+'.png')
  plt.close()

# Plot ice timesnaps
for i,t in enumerate(time_ice):

  if i % 20 == 0:
    plt.figure()
    plt.contourf(lon_ice,lat_ice,ice_nc.variables['ICEC_L1'][i,:,:])
    plt.savefig('ice_'+str(i).zfill(4)+'.png')
    plt.close()

