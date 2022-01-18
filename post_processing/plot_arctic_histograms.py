import matplotlib.pyplot as plt
import numpy as np
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from geometric_features import read_feature_collection
import shapely.geometry
from matplotlib import cm
import pprint
import re

regions = [
           'Chukchi_Sea_NSIDC',
           'Bering_Sea',
           'Beaufort_Sea_NSIDC',
           'East_Siberian_Sea_NSIDC',
           'Canadian_Archipelago_NSIDC',
           'Laptev_Sea_NSIDC',
           'Hudson_Bay_NSIDC',
           'Kara_Sea',
           'Baffin_Bay_NSIDC',
           'Barents_Sea',
           'Greenland_Sea',
           'Central_Arctic_NSIDC',
          ]
regions_path = '/Users/sbrus/Data/geometric_features/geometric_data/ocean/region/'


seasons = ['ANN','JFM','JAS']

fields = ['significantWaveHeight','peakWaveFrequency']

ylimits = {'significantWaveHeight':{'field':[0.0,1.0],'diff':[-0.2,0.2]},
           'peakWaveFrequency':{'field':[0.0,0.4],'diff':[-0.1,0.1]}}

units = ['(m)','(s)']

year = 2090
year_control = 2020
base_path = '/Users/sbrus/Projects/E3SM/20220106.AWAV_WCYCLSSP585_CMIP6.ne30_ECv3_wQU225EC60to30.chrysalis/analysis_'
data_path = base_path+str(year)+'/tables/waves/'
data_path_control = base_path+str(year_control)+'/tables/waves/'

#lightness
#cmap = ["#5a3a8e",
#"#b0457b",
#"#b54e36",
#"#ba4758",
#"#6d8838",
#"#6d83da",
#"#c771c4",
#"#ce822e",
#"#bb8c45",
#"#5fbb69",
#"#a3b23f",
#"#45c097"]

cmap = ["#5a3a8e",
"#a3b23f",
"#b0457b",
"#6d8838",
"#6d83da",
"#45c097",
"#ba4758",
"#5fbb69",
"#c771c4",
"#bb8c45",
"#b54e36",
"#ce822e"]

offset = {
           'Chukchi_Sea_NSIDC':[0.5,1.0],
           'Bering_Sea':[7.0,0.0],
           'Beaufort_Sea_NSIDC':[2.0,2.0],
           'East_Siberian_Sea_NSIDC':[0.0,2.0],
           'Canadian_Archipelago_NSIDC':[-5.0,-10.0],
           'Laptev_Sea_NSIDC':[-2.0,2.5],
           'Hudson_Bay_NSIDC':[-1.0,-1.0],
           'Kara_Sea':[3.0,2.0],
           'Baffin_Bay_NSIDC':[7.0,0.0],
           'Barents_Sea':[-3.0,0.0],
           'Greenland_Sea':[5.0,0.0],
           'Central_Arctic_NSIDC':[0.0,5.0],
          }

########################################################################
########################################################################

def read_hist_data(path):
  f = open(path,'r')
  x = []
  w = 0
  h = []
  lines = f.readlines()
  for line in lines:
    line_sp = line.split()
    x.append(float(line_sp[0]))
    w = float(line_sp[1])
    h.append(float(line_sp[2]))

  x = np.array(x)
  h = np.array(h)

  return x,w,h

def plot_regions(ax):

  ax.set_extent([0, 359.9, 50, 90], crs=ccrs.PlateCarree())
  ax.add_feature(cfeature.LAND, alpha=1, zorder=100)
  #ax.add_feature(cfeature.LAKES, alpha=1, zorder=101)
  ax.add_feature(cfeature.COASTLINE, zorder=101)
  ax.gridlines(linestyle='dotted')
  
  for i,region in enumerate(regions):
  
    path = regions_path+region 
    fc = read_feature_collection(path+'/region.geojson')
    #pprint.pprint(fc.features)
  
    shape = []
    for s in fc.features[0]['geometry']['coordinates']:
      if len(s) == 1:
        shp = shapely.geometry.Polygon(s[0])
      else:
        shp = shapely.geometry.Polygon(s)
      shape.append(shp)
    
    ax.add_geometries(shape,crs=ccrs.PlateCarree(),facecolor=cmap[i])
    center = list(shape[0].centroid.coords)[0]
  
    ax.text(center[0]+offset[region][0],
             center[1]+offset[region][1],
             region.replace('_NSIDC','').replace('_','\n'),
             transform=ccrs.PlateCarree(),
             ha='center',va='center',zorder=102,
             fontsize=9)

########################################################################
########################################################################

fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.NorthPolarStereo())
plot_regions(ax)
fig.tight_layout()
fig.savefig('regions.png')
plt.close()

for diff in [True,False]:

  for j,field in enumerate(fields):
    print(field)
    for season in seasons:
      print('  '+season)
  
      fig2,axd = plt.subplot_mosaic([['.','Chukchi_Sea_NSIDC','Bering_Sea','.'],
                                    ['Beaufort_Sea_NSIDC',        'region','region','East_Siberian_Sea_NSIDC'],
                                    ['Canadian_Archipelago_NSIDC','region','region','Laptev_Sea_NSIDC'],
                                    ['Hudson_Bay_NSIDC',          'region','region','Kara_Sea'],
                                    ['Baffin_Bay_NSIDC',          'region','region','Barents_Sea'],
                                    ['.',                         'Greenland_Sea','Central_Arctic_NSIDC','.']],figsize=(12,9))

      xx = axd['region'].get_subplotspec()
      axd['region'].remove()
      axd['region'] = fig2.add_subplot(xx,projection=ccrs.NorthPolarStereo())
      plot_regions(axd['region'])

      for i,region in enumerate(regions):
        print('    '+region)
        
        filename = field + '_' + season + '_'+ region + '.txt'
  
        path = data_path + filename
        x,w,h = read_hist_data(path)
        ylim = ylimits[field]['field']

        year_str = str(year)
  
        path = data_path_control + filename
        x,w,h_control = read_hist_data(path)
        if diff:
          h = h-h_control
          year_str = year_str+'_control'+str(year_control)
          ylim = ylimits[field]['diff'] 
          
  
        axd[region].bar(x,h,width=w,align='edge',color=cmap[i],edgecolor='k')
        if not diff:
          axd[region].bar(x,h_control,width=w,align='edge',alpha=0.6,color='grey',edgecolor='k')
        axd[region].set_ylim(ylim)
        region_name = region.replace('_NSIDC','').replace('_',' ')
        axd[region].set_title(region_name,fontsize=9)
  
      field_name_sp = re.split('(?=[A-Z])',field)
      field_name_sp[0] = field_name_sp[0].capitalize()
      field_name = ' '.join(field_name_sp)+' '+units[j]

      ax = fig2.add_subplot(111,frameon=False)
      ax.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
      ax.set_xlabel(field_name,fontsize=14,labelpad=5)
      ax.set_ylabel('Density',fontsize=14,labelpad=10)

      fig2.tight_layout()
      fig2.savefig(field+'_'+season+'_'+year_str+'_All_Regions.png')
      plt.close()

