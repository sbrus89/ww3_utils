import numpy as np
import xarray as xr
import numpy.ma as ma
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
import matplotlib.ticker as ticker
import netCDF4 
import glob
import datetime
import calendar
import os
import yaml
import pprint
import subprocess
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy import interpolate
import pprint
plt.switch_backend('agg')


###############################################################################################
###############################################################################################

def read_field_ww3(f,sol,variable='hs'):

  nc_file = netCDF4.Dataset(f,'r')

  # Handle date
  ref_date = nc_file.variables['time'].getncattr('units').replace('days since ','')
  ref_date = datetime.datetime.strptime(ref_date,'%Y-%m-%d %H:%M:%S')
  output_time = nc_file.variables['time'][:]
  if len(output_time) > 1:
    print("Should be only one timestep per nc file. Check ww3_ounf.inp")
    raise SystemExit(0)
  sol['datetime'] = ref_date + datetime.timedelta(days=output_time[0])
  sol['date'] = sol['datetime'].strftime('%Y %m %d %H %M')

  # Mesh type  
  if 'element' in nc_file.dimensions:
    sol['type'] = 'unstruc'
  else:
    sol['type'] = 'struc'

  # Decide if mesh should be read
  if 'lon' not in sol:
    read_mesh = True
  else:
    read_mesh = False

  # Read mesh, if necessary 
  if read_mesh:
    sol['lon']  = nc_file.variables['longitude'][:] 
    sol['lat']  = nc_file.variables['latitude'][:]
    sol['mapsta'] = nc_file.variables['MAPSTA'][:]
    if sol['type'] == 'struc':
      idx = np.where(sol['lon'] > 180.0)
      sol['lon'][idx] = sol['lon'][idx] - 360.0
      idx = np.argsort(sol['lon'])
      sol['lon'] = sol['lon'][idx]
      sol['idx'] = idx
  if read_mesh and sol['type'] == 'unstruc':
    sol['tri'] = nc_file.variables['tri'][:] - 1
    sol['tri'] = Triangulation(sol['lon'],sol['lat'],sol['tri'])
   
  # Read variable
  if sol['type'] == 'struc':
    sol[variable] = nc_file.variables[variable][0,:,:]
    sol[variable] = sol[variable][:,sol['idx']]
  else:
    sol[variable] = nc_file.variables[variable][0,:]

  nc_file.close()

###############################################################################################
###############################################################################################

def interpolate_solution(sol_interp,sol_grid,variable='hs'):
 
  if sol_interp['type'] == 'unstruc': 
    interp = interpolate.LinearNDInterpolator((sol_interp['lon'],sol_interp['lat']),sol_interp[variable])
  else:
    interp = interpolate.RegularGridInterpolator((sol_interp['lon'],sol_interp['lat']),sol_interp[variable].T, bounds_error=False)

  Lon,Lat = np.meshgrid(sol_grid['lon'],sol_grid['lat'])
  pts = np.vstack((Lon.ravel(),Lat.ravel())).T
  sol_interp[variable+'_interp'] = interp(pts)
  sol_interp[variable+'_interp'] = np.reshape(sol_interp[variable+'_interp'],Lon.shape)
  sol_interp[variable+'_interp'] = ma.masked_array(sol_interp[variable+'_interp'],mask=sol_grid[variable].mask)

###############################################################################################
###############################################################################################

def plot_differences(sol_ref,sol_interp,variable='hs'):

  levels = np.linspace(0,10,100)
  fig = plt.figure(figsize=(6,14))

  # Uninterpolated solution
  ax = fig.add_subplot(5, 1, 1,projection=ccrs.PlateCarree())
  if sol_interp['type'] == 'unstruc':
    sol_orig = ma.filled(sol_interp[variable],fill_value=0.0)
    cf = ax.tricontourf(sol_interp['lon'],sol_interp['lat'],sol_orig, levels=levels, transform=ccrs.PlateCarree())
  else:
    cf = ax.contourf(sol_interp['lon'],sol_interp['lat'],sol_interp[variable], levels=levels, transform=ccrs.PlateCarree())
  ax.add_feature(cfeature.LAND, zorder=100)
  ax.add_feature(cfeature.LAKES, alpha=0.5, zorder=101)
  ax.add_feature(cfeature.COASTLINE, zorder=101)
  plt.colorbar(cf,ax=ax)

  # Interpolated solution onto reference solution grid
  ax = fig.add_subplot(5, 1, 2,projection=ccrs.PlateCarree())
  cf = ax.contourf(sol_ref['lon'],sol_ref['lat'],sol_interp[variable+'_interp'], levels=levels, transform=ccrs.PlateCarree())
  ax.add_feature(cfeature.LAND, zorder=100)
  ax.add_feature(cfeature.LAKES, alpha=0.5, zorder=101)
  ax.add_feature(cfeature.COASTLINE, zorder=101)
  plt.colorbar(cf,ax=ax)

  # Reference solution
  ax = fig.add_subplot(5, 1, 3,projection=ccrs.PlateCarree())
  cf = ax.contourf(sol_ref['lon'],sol_ref['lat'],sol_ref[variable], levels=levels, transform=ccrs.PlateCarree())
  ax.add_feature(cfeature.LAND, zorder=100)
  ax.add_feature(cfeature.LAKES, alpha=0.5, zorder=101)
  ax.add_feature(cfeature.COASTLINE, zorder=101)
  plt.colorbar(cf,ax=ax)

  # Absolute Difference
  levels = np.linspace(-1.0,1.0,100)
  ax = fig.add_subplot(5, 1, 4,projection=ccrs.PlateCarree())
  cf = ax.contourf(sol_struc['lon'],sol_struc['lat'],sol_interp[variable+'_interp']-sol_struc[variable], levels=levels, cmap='bwr', transform=ccrs.PlateCarree())
  ax.add_feature(cfeature.LAND, zorder=100)
  ax.add_feature(cfeature.LAKES, alpha=0.5, zorder=101)
  ax.add_feature(cfeature.COASTLINE, zorder=101)
  plt.colorbar(cf,ax=ax)

  # Relative Difference
  levels = np.linspace(-100.0,100.0,100)
  ax = fig.add_subplot(5, 1, 5,projection=ccrs.PlateCarree())
  cf = ax.contourf(sol_struc['lon'],sol_struc['lat'],100.0*np.divide(sol_interp[variable+'_interp']-sol_struc[variable],sol_struc[variable]), levels=levels, cmap='bwr', transform=ccrs.PlateCarree())
  ax.add_feature(cfeature.LAND, zorder=100)
  ax.add_feature(cfeature.LAKES, alpha=0.5, zorder=101)
  ax.add_feature(cfeature.COASTLINE, zorder=101)
  plt.colorbar(cf,ax=ax)

  plt.savefig('interp.png', bbox_inches='tight', dpi=400)
  plt.close()

###############################################################################################
###############################################################################################

def plot_stats(lon,lat,mean,maximum,title,filename):

  fig = plt.figure(figsize=(6,6))

  # Mean Difference
  levels = np.linspace(0.0,0.5,100)
  ax = fig.add_subplot(2, 1, 1,projection=ccrs.PlateCarree())
  cf = ax.contourf(lon,lat,mean, levels=levels, cmap='RdYlGn_r', extend='max',transform=ccrs.PlateCarree())
  ax.add_feature(cfeature.LAND, zorder=100)
  ax.add_feature(cfeature.LAKES, alpha=0.5, zorder=101)
  ax.add_feature(cfeature.COASTLINE, zorder=101)
  ax.set_title('Mean Difference\n'+title)
  cb = plt.colorbar(cf,ax=ax,ticks=np.linspace(0.0,0.5,6),format=ticker.FormatStrFormatter('%4.2f'))
  cb.set_label('mean difference (m)')
  
  # Max Difference
  levels = np.linspace(0.0,1.0,100)
  ax = fig.add_subplot(2, 1, 2,projection=ccrs.PlateCarree())
  cf = ax.contourf(lon,lat,maximum, levels=levels, cmap='RdYlGn_r', extend='max', transform=ccrs.PlateCarree())
  ax.add_feature(cfeature.LAND, zorder=100)
  ax.add_feature(cfeature.LAKES, alpha=0.5, zorder=101)
  ax.add_feature(cfeature.COASTLINE, zorder=101)
  ax.set_title('Maximum Difference\n'+title)
  cb = plt.colorbar(cf,ax=ax,ticks=np.linspace(0.0,1.0,6),format=ticker.FormatStrFormatter('%4.2f'))
  cb.set_label('maximum difference (m)')
  

  plt.savefig('diff_stats_'+filename+'.png', bbox_inches='tight', dpi=400)
  plt.close()

###############################################################################################
###############################################################################################

if __name__ == '__main__':

  #run_unst = 'unst_depth_4000_UOST_edit2'
  #run_unst = 'glo_2d_UOST_PR1'
  run_unst = 'unst_depth_4000'
  run_struc = 'glo_30m_UOST_PR1'

  dir_unst = '/users/sbrus/scratch4/WW3_testing/'+run_unst+'/post_processing/june/model_data/fields/'
  dir_struc = '/users/sbrus/scratch4/WW3_testing/'+run_struc+'/post_processing/june/model_data/fields/'

  if run_unst == 'unst_depth_4000_UOST_edit2':
    title = '(unstructured - 1/2 degree structured)'
  if run_unst == 'unst_depth_4000':
    title = '(unstructured no UOST - 1/2 degree structured)'
  elif run_unst == 'glo_2d_UOST_PR1':
    title = '(2 degree structured - 1/2 degree structured)'
  
  file_paths = sorted(glob.glob(dir_unst+'*.nc'))
  files = []
  for f in file_paths:
    fname = f.split('/')[-1]
    files.append(fname)
  

  for i,f in enumerate(files):
    print(f)

    sol_unst = {}
    read_field_ww3(dir_unst+f,sol_unst)
    #pprint.pprint(sol_unst)

    sol_struc = {}
    read_field_ww3(dir_struc+f,sol_struc)
    #pprint.pprint(sol_struc)

    interpolate_solution(sol_unst,sol_struc)

    if i == 0:
      mean = np.zeros(sol_struc['hs'].shape)
      maximum = np.zeros(sol_struc['hs'].shape)

    d = sol_unst['hs_interp']-sol_struc['hs']
    mean = mean + np.absolute(d)
    maximum = np.maximum(maximum,np.absolute(d))

  mean = mean/float(len(files))

  filename = run_unst+'-'+run_struc
  plot_stats(sol_struc['lon'],sol_struc['lat'],mean,maximum,title,filename)


  #plot_differences(sol_struc,sol_unst)
