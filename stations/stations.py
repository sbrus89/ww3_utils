import xml.etree.ElementTree as ET
import os
import urllib2
import datetime
import pprint

pwd = os.getcwd()

# Input parameters
year = '2005'
run_start_date = '01-01'
run_end_date   = '12-31'
movement_tolerance = 0.1

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
    ID = "'"+sta.attrib['id']+"'"
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

      # Determine if station has data in the date range for the run
      if (run_start <= sta_end and run_end > sta_start) and (met == 'y'):
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
        

n = 0
for sta in stations:
  if stations[sta]['moved'] == True:
    n = n + 1
    pprint.pprint(stations[sta])
print n    

#pprint.pprint(stations)
#print len(stations)


f = open(pwd+'/stations.txt','w')
for sta in stations:
 

      if stations[sta]['moved'] == True:
        continue

      owner = stations[sta]['owner']
      prgm = stations[sta]['prgm']
      lon = stations[sta]['history'][0]['lon']
      lat = stations[sta]['history'][0]['lat']

      success = False
      try:
        print 'downloading '+sta
        url = 'https://www.ndbc.noaa.gov/view_text_file.php?fileID='+sta+'h'+year+'.txt.gz&dir=data/historical/stdmet/'
        data = urllib2.urlopen(url).read().splitlines()
        success = True
      except:
        print '  error downloading data'

      if success:
        print data[0]
        if data[0].find('YYYY MM DD hh mm') >= 0:

          # Save station data file
          fd = open(sta.attrib['ID']+'_'+year+'.txt','w')
          fd.write('\n'.join(data))
          fd.close()

          # Add to stations list file
          f.write('  '.join([lon,lat,sta,owner,prgm])+'\n')
        else: 
          print '  data not availiable'
    
