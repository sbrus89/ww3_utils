import xml.etree.ElementTree as ET
import os
import urllib2
import datetime
import pprint
import glob

# Creates station file from a set of downloaded data files


pwd = os.getcwd()

# Input parameters
year = '2009'
run_start_date = '01-01'
run_end_date   = '12-31'
movement_tolerance = 0.1
lonlat_box = [-180.0,180.0,-90.0,90.0]
data_product = 'standard meterological'
#data_product = 'spectral wave'

data_product_ID = {'standard meterological':[['stdmet','h']],
                   'spectral wave':[['swden','w'],['swdir','d'],['swdir2','i'],['swr1','j'],['swr2','k']]}

# Convert dates to datetime objects
datetime_frmt = '%Y-%m-%d'
run_start = datetime.datetime.strptime(year+'-'+run_start_date,datetime_frmt)
run_end   = datetime.datetime.strptime(year+'-'+run_end_date,datetime_frmt)

# Open and parse station data file
fdir = os.path.dirname(os.path.realpath(__file__))
tree = ET.parse(fdir+'/stationmetadata.xml')
root = tree.getroot()

stations = {}
nsta = 0
for i,sta in enumerate(root):

    # Get station attribute data
    ID = "'"+sta.attrib['id'].lower()+"'"
    owner = sta.attrib['owner']
    prgm = sta.attrib['pgm']      
    #print i+1, ID, prgm

    for hist in sta:

      # Get station history data
      lon = hist.attrib['lng']
      lat = hist.attrib['lat']
      met = hist.attrib['met']
      sta_start_date = hist.attrib['start']
      sta_end_date = hist.attrib['stop']

      # Convert dates to datetime objects
      if sta_end_date == '':
        sta_end_date = datetime.datetime.today().strftime(datetime_frmt)
      sta_start = datetime.datetime.strptime(sta_start_date,datetime_frmt)
      sta_end = datetime.datetime.strptime(sta_end_date,datetime_frmt)

      # Determine if stations is inside lon/lat box:w
      inside_box = False
      if float(lon) >= lonlat_box[0] and float(lon) <= lonlat_box[1] and float(lat) >= lonlat_box[2] and float(lat) <= lonlat_box[3]:
        inside_box = True

      # Determine if station has data in the date range for the run
      if (run_start <= sta_end and run_end > sta_start) and (met == 'y') and inside_box:
        if ID not in stations:
          stations[ID] = {}
          stations[ID]['owner'] = owner
          stations[ID]['prgm'] = prgm
          stations[ID]['history'] = []
          stations[ID]['moved'] = False
          nsta = nsta + 1
          print nsta,ID,prgm
        stations[ID]['history'].append({'start':sta_start_date,'end':sta_end_date,'lon':lon,'lat':lat})
        print '  ', hist.attrib['start']

    # Check if station location has moved throughout deployment
    if (ID in stations) and (len(stations[ID]['history']) > 1):
      moved = False
      for hist in stations[ID]['history']:
        lon_movement = abs(float(hist['lon'])-float(stations[ID]['history'][0]['lon']))
        lat_movement = abs(float(hist['lat'])-float(stations[ID]['history'][0]['lat']))
        if lon_movement > movement_tolerance or \
           lat_movement > movement_tolerance:
          moved = True
      if moved == True:
        stations[ID]['moved'] = True

        
files = glob.glob(pwd+'/*_'+year+'_*.txt')
print files

# Download station data and write to station list file
f = open(pwd+'/stations_test.txt','w')

sta_list = []
for sta_file in files: 
  
  sta = "'"+sta_file.split('/')[-1].split('_')[0]+"'"
  print sta

  # Get station information
  owner = stations[sta]['owner']
  prgm = stations[sta]['prgm']
  lon = stations[sta]['history'][0]['lon']
  lat = stations[sta]['history'][0]['lat']

  # Add to stations list file
  if sta not in sta_list:
    f.write('  '.join([lon,lat,sta,owner,prgm])+'\n')
    sta_list.append(sta) 

f.close()
