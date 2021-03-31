import os
import numpy as np
from scipy import interpolate
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import calc_distance
plt.switch_backend('agg')
cartopy.config['pre_existing_data_dir'] = \
    os.getenv('CARTOPY_DIR', cartopy.config.get('pre_existing_data_dir'))

km = 1000.0

depth_threshold_refined = 1000.0
distance_threshold_refined = 300.0
depth_threshold_global= 1000.0
distance_threshold_global= 100.0
refined_res = 20*km

maxres = 225.0*km

def cell_widthVsLatLon(lon,lat,coast_shapefiles,sphere_radius,ocean_mesh):


    # Create structrued grid points
    nlon = lon.size
    nlat = lat.size
    Lon_grd, Lat_grd = np.meshgrid(np.radians(lon), np.radians(lat))
    xy_pts = np.vstack((Lon_grd.ravel(),Lat_grd.ravel())).T

    # Create structured grid points and background cell_with array 
    cell_width = np.zeros(Lon_grd.shape) + maxres

    # Get ocean mesh variables
    nc_file = Dataset(ocean_mesh,'r')
    areaCell = nc_file.variables['areaCell'][:]
    lonCell = nc_file.variables['lonCell'][:]
    latCell = nc_file.variables['latCell'][:]
    bottomDepth = nc_file.variables['bottomDepth'][:]

    # Transform 0,360 range to -180,180 
    idx, = np.where(lonCell > np.pi) 
    lonCell[idx] = lonCell[idx] - 2.0*np.pi
    idx, = np.where(lonCell < -np.pi)
    lonCell[idx] = lonCell[idx] + 2.0*np.pi

    # Interpolate cellWidth onto background grid
    cellWidth = 2.0*np.sqrt(areaCell/np.pi)
    hfun = interpolate.NearestNDInterpolator((lonCell,latCell),cellWidth)
    hfun_interp = hfun(xy_pts)
    hfun_grd = np.reshape(hfun_interp,(nlat,nlon))

    # Interpolate bathymetry onto background grid
    bathy = interpolate.NearestNDInterpolator((lonCell,latCell),bottomDepth)
    bathy_interp = bathy(xy_pts)
    bathy_grd = np.reshape(bathy_interp,(nlat,nlon))

    # Get distance to coasts
    D = calc_distance.distance_to_shapefile_points(coast_shapefiles,lon,lat,sphere_radius,reggrid=True)

    # Apply refined region criteria
    idxx,idxy = np.where(   (bathy_grd < depth_threshold_refined) 
                          & (bathy_grd > 0.0) 
                          & (D < distance_threshold_refined) 
                          & (hfun_grd < refined_res))
    cell_width[idxx,idxy] = hfun_grd[idxx,idxy]
    
    # Apply global region criteria
    idxx,idxy = np.where(   (bathy_grd < depth_threshold_global) 
                          & (bathy_grd > 0.0) 
                          & (D < distance_threshold_global))
    cell_width[idxx,idxy] = hfun_grd[idxx,idxy]


    # Plot  
    fig = plt.figure(figsize=(16,8))

    ax1 = fig.add_subplot(4, 1, 1, projection=ccrs.PlateCarree())
    plt1 = ax1.contourf(lon, lat, cell_width, transform=ccrs.PlateCarree())
    ax1.set_title('waves mesh cell width')
    fig.colorbar(plt1,ax=ax1)

    ax2 = fig.add_subplot(4, 1, 2, projection=ccrs.PlateCarree())
    plt2 = ax2.contourf(lon, lat, hfun_grd, transform=ccrs.PlateCarree())
    ax2.set_title('ocean cell width')
    fig.colorbar(plt2,ax=ax2)

    ax3 = fig.add_subplot(4, 1, 3, projection=ccrs.PlateCarree())
    plt3 = ax3.contourf(lon, lat, bathy_grd, transform=ccrs.PlateCarree())
    ax3.set_title('bathymetry')
    fig.colorbar(plt3,ax=ax3)

    ax4 = fig.add_subplot(4, 1, 4, projection=ccrs.PlateCarree())
    plt4 = ax4.contourf(lon, lat, D, transform=ccrs.PlateCarree())
    ax4.set_title('distance to coast')
    fig.colorbar(plt4,ax=ax4)

    fig.savefig(
        'cell_width' +
        '.png',
        bbox_inches='tight',
        dpi=400)
    plt.close()

    return cell_width 

