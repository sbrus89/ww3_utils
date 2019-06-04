import xml.etree.ElementTree as ET
import os
import urllib2
import datetime
import pprint

pwd = os.getcwd()

# Input parameters
year = '1999'
run_start_date = '01-01'
run_end_date   = '12-31'
movement_tolerance = 0.1
lonlat_box = [-180.0,180.0,-90.0,90.0]
data_product = 'standard meterological'

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

        
## Print stations that have moved
#n = 0
#for sta in stations:
#  if stations[sta]['moved'] == True:
#    n = n + 1
#    pprint.pprint(stations[sta])
#print n    

#pprint.pprint(stations)
#print len(stations)


#Ignore proxies (for use outside lab)
#urllib2.getproxies = lambda: {}


# Download station data and write to station list file
f = open(pwd+'/stations.txt','w')
for sta in stations:
      
      print "Station "+sta
 
      # Skip if station has moved
      if stations[sta]['moved'] == True:
        print " station has moved"
        continue

      # Get station information
      owner = stations[sta]['owner']
      prgm = stations[sta]['prgm']
      lon = stations[sta]['history'][0]['lon']
      lat = stations[sta]['history'][0]['lat']
      ID = sta.strip("'")

      # Additional check to see if all data products exist
      all_data_exists = True
      station_found = True
      for prod in data_product_ID[data_product]:
        try:
          print '  checking '+ID+' '+prod[0]
          url = 'https://www.ndbc.noaa.gov/station_history.php?station='+ID
          data = urllib2.urlopen(url).read()
          if data.find(ID+prod[1]+year) < 0:
            all_data_exists = False
            break
        except:
          station_found = False

      if not station_found:
        print '  station not found'
        continue 

      if not all_data_exists:
        print "  not all data products exist for this year"
        continue
 
      # Download data
      success = False
      added_to_list = False
      for prod in data_product_ID[data_product]:
        try:
          print '  downloading '+ID+ ' '+prod[0]

          url = 'https://www.ndbc.noaa.gov/view_text_file.php?filename='+ID+prod[1]+year+'.txt.gz&dir=data/historical/'+prod[0]+'/'
          data = urllib2.urlopen(url).read().splitlines()
          success = True
        except:
          print '  error downloading data'

        # Write to file
        if success:
          print data[0]
          if data[0].find('YYYY MM DD hh mm') >= 0 or data[0].find('#YY  MM DD hh mm') >= 0 or data[0].find('YYYY MM DD hh') >= 0 or data[0].find('YY MM DD hh') >= 0:

            # Save station data file
            fd = open(ID+'_'+year+'_'+prod[0]+'.txt','w')
            fd.write('\n'.join(data))
            fd.close()

            # Add to stations list file
            if not added_to_list:
              f.write('  '.join([lon,lat,sta,owner,prgm])+'\n')
              added_to_list = True
          else: 
            print '  data not availiable'
    
