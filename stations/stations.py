import xml.etree.ElementTree as ET
import os
import urllib2


pwd = os.getcwd()
fdir = os.path.dirname(os.path.realpath(__file__))


tree = ET.parse(fdir+'/activestations.xml')
root = tree.getroot()

owners = ['NDBC']
year = '2005'

f = open(pwd+'/stations.txt','w')
for sta in root:
  for owner in owners:
    if sta.attrib['owner'] == owner:
      lon = sta.attrib['lon']
      lat = sta.attrib['lat']
      name = "'"+sta.attrib['id']+"'"
      prgm = sta.attrib['pgm']      
      f.write('  '.join([lon,lat,name,owner,prgm])+'\n')
      
      success = False
      try:
        print 'downloading '+sta.attrib['id']
        url = 'https://www.ndbc.noaa.gov/view_text_file.php?filename='+sta.attrib['id']+'h'+year+'.txt.gz&dir=data/historical/stdmet/'
        data = urllib2.urlopen(url).read().splitlines()
        success = True
      except:
        print '  error downloading data'

      if success:
        print data[0]
        if data[0].find('YYYY MM DD hh mm') >= 0:
          fd = open(sta.attrib['id']+'_'+year+'.txt','w')
          fd.write('\n'.join(data))
          fd.close()
        else: 
          print '  data not availiable'
     
