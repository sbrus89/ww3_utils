import xml.etree.ElementTree as ET
import os

pwd = os.getcwd()
fdir = os.path.dirname(os.path.realpath(__file__))


tree = ET.parse(fdir+'/activestations.xml')
root = tree.getroot()

owners = ['NDBC']

f = open(pwd+'/stations.txt','w')
for sta in root:
  for owner in owners:
    if sta.attrib['owner'] == owner:
      lon = sta.attrib['lon']
      lat = sta.attrib['lat']
      name = "'"+sta.attrib['id']+"'"
      prgm = sta.attrib['pgm']      
      f.write('  '.join([lon,lat,name,owner,prgm])+'\n')



