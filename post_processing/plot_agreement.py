import plot_points
from scipy import stats
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.cm as cm
import matplotlib.colors as colors
import glob
import datetime
import numpy as np
import os

direc = '/users/sbrus/scratch4/WW3_testing/glo_30m/post_processing/june-august/model_data/points/'
obs_direc = './obs_data/'
wav_files = [sorted(glob.glob(direc+'ww3*_tab.nc'))]
variables = plot_points.variables
variable = 'hs'
year = '2005'

# Read in modeled data
data, stations, output_time, ref_date = plot_points.read_point_files(wav_files,variables)

# Convert output times to date format
output_date = []
output_datetime = []
for t in output_time:
  date = ref_date + datetime.timedelta(days=t)
  output_date.append(date.strftime('%Y %m %d %H %M'))
  output_datetime.append(date)
output_datetime = np.asarray(output_datetime,dtype='O')

# Initialize plot
fig = plt.figure(figsize=[12,12])
m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
                             llcrnrlon=-180,urcrnrlon=-30,resolution='c')
m.fillcontinents(color='tan',lake_color='lightblue')
m.drawcoastlines()
#m.drawparallels(np.arange(-90,90,20),labels=[False,True,True,False])
#m.drawmeridians(np.arange(-180,180,20),labels=[True,False,False,True])
cmap = cm.RdYlGn_r
norm = colors.Normalize(vmin=0.0,vmax=100)


xv = []
yv = []
zv = []
for sta in sorted(stations['name']):
  print sta

  # Get data from observation file at output times
  obs_file = obs_direc+sta+'_'+year+'.txt'
  if os.path.isfile(obs_file):
    obs_data = plot_points.read_station_data(obs_file,output_date,variables)

    ind = stations['name'].index(sta)
    lon = stations['lon'][ind]
    lat = stations['lat'][ind]

    # Only plot stations with modeled data (wet points)
    if len(np.unique(data[variable][:,ind])) > 1:

      mask = np.isnan(obs_data[variable])
      idx, = np.where(mask == False)
      x = data[variable][:,ind][idx]
      y = obs_data[variable][idx]
      
      # Only plot stations with observed data
      if len(idx) > 1:

        # Calculate regression coeffieint
#        slope,b,r,p,err = stats.linregress(x,y)
#        r2 = r**2
#        print r2

        r2 = np.mean(np.absolute(np.divide(y-x,y)))*100.0

        xv.append(lon)
        yv.append(lat)
        zv.append(r2)

colors = [cmap(norm(val)) for val in zv]
sc = m.scatter(xv,yv,c=zv,marker='o',cmap=cmap,zorder=10,vmin=0.0,vmax=100.0)
cbar = plt.colorbar(sc,orientation='horizontal')
#cbar.set_label('r squared value')
cbar.set_label('Mean absolute percentage error')
fig.tight_layout()
fig.savefig('map.png',bbox_inches='tight',dpi=400)    
