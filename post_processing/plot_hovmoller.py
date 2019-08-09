import glob
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from plot_fields import read_field_ww3 
from geopy import distance
from scipy import interpolate
import datetime
plt.switch_backend('agg')

if __name__ == '__main__':
  
  ww3_direcs = ['/users/sbrus/scratch4/WW3_CFSR_2000-2010/2002/run00/results/model_data/fields/']
  transects = [
               {'lon0':-65.0, 'lat0':-65.0, 'lon1':-65.0, 'lat1':-55.0, 'spacing':10.0, 'name':'south america'},
               {'lon0': 20.0, 'lat0':-69.0, 'lon1': 20.0, 'lat1':-35.0, 'spacing':10.0, 'name':'africa'},
               {'lon0':129.0, 'lat0':-66.0, 'lon1':129.0, 'lat1':-31.5, 'spacing':10.0, 'name':'austrailia'}
              ]

  # Create list of ww3 field output files
  ww3_files = []
  for direc in ww3_direcs:  
    ww3_files.extend(glob.glob(direc+'*.nc'))
  ww3_files.sort()
  nfiles = len(ww3_files)

  # Create list of transect data
  transect_points = []
  for xsec in transects:
    tot_dist = distance.distance((xsec['lat0'],xsec['lon0']),(xsec['lat1'],xsec['lon1'])).km
    print tot_dist
    n = int(tot_dist/xsec['spacing'])
    dist = np.linspace(0,tot_dist,n)
    lon = np.zeros((1,n))
    lon = np.linspace(xsec['lon0'],xsec['lon1'],n)
    lat = np.zeros((1,n))
    lat = np.linspace(xsec['lat0'],xsec['lat1'],n)
    data = np.zeros((nfiles,n))
    dates = []
    transect_points.append({'lon':lon,'lat':lat,'data':data,'dates':dates,'dist':dist,'name':xsec['name']})

  
  for n,f in enumerate(ww3_files):
    print f
    lon_vec,lat_vec,hs,output_date = read_field_ww3(f)

    for xsec in transect_points:
      interp = interpolate.RegularGridInterpolator((lon_vec, lat_vec), hs.T, bounds_error=False, fill_value=np.nan)    
      lon = np.mod(xsec['lon']+360,360)
      pts = np.vstack((lon,xsec['lat'])).T
      xsec['data'][n,:] = interp(pts)
      xsec['dates'].append(datetime.datetime.strptime(output_date,'%Y %m %d %H %M'))


  for xsec in transect_points:

    fig = plt.figure(figsize=[12.0,4.0])
    ax = fig.add_subplot(2,1,1)
    m = Basemap(projection='cyl',llcrnrlat= -90,urcrnrlat=90,\
                                 llcrnrlon=-180,urcrnrlon=180,resolution='c')
    m.fillcontinents(color='tan',lake_color='lightblue')
    m.drawcoastlines()
    ax.plot(xsec['lon'],xsec['lat'],'r-')

    ax = fig.add_subplot(2,1,2)
    dates = np.asarray(xsec['dates'],dtype='O')
    print xsec['data'].shape
    print dates.shape
    print xsec['dist'].shape
    idx = np.where(np.absolute(xsec['data']) > 1e10)
    xsec['data'][idx] = np.nan
    cf = ax.contourf(dates,xsec['dist'],xsec['data'].T)
    cb = plt.colorbar(cf)
    cb.set_label('Significant Wave Height')
    plt.tight_layout()
    plt.savefig('hovmoller_'+'_'.join(xsec['name'].split())+'.png',bbox_inches='tight')
      
