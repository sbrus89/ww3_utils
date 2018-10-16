import subprocess
import os
import glob

# Combines two CFSR NetCDF files for a given month
# into a single NetCDF file that can be read by ww3_prnc.
# 
# For whatever reason, the hourly wind data for a given month
# is split into two files.
#
# The first  file contains only the 00 hour data for the first 
#    day of the month contained in the second file.
# The second file contains the rest of the hourly data for that month, 
#    beginning with the 01 hour of the first day.
#
# The latitude dimension also has to be altered for ww3_prnc
#
# For ocean current data, the u and v components must first be
# appended into a single file.



pwd = os.getcwd()+'/'

files = sorted(glob.glob(pwd+'*.nc'))

first_file  = files[0].split('/')[-1] 
second_file = files[1].split('/')[-1]

if files[0].find('ocn') >= 0:
  print 'Appending files...'  
  fileA = files[0].split('/')[-1]
  fileB = files[2].split('/')[-1]
  subprocess.call(['ncks -A',pwd+fileB,pwd+fileA])
  fileA = files[1].split('/')[-1]
  fileB = files[3].split('/')[-1]
  subprocess.call(['ncks -A',pwd+fileB,pwd+fileA])

print first_file
print second_file

print 'Concatenatng files...'
subprocess.call(['ncrcat',pwd+first_file,second_file,pwd+'tmp.nc'])
print 'Fixing latitude dimension...'
subprocess.call(['ncpdq','-O','-a','-lat',pwd+'tmp.nc',pwd+second_file.replace('grb2','ww3')])
print 'Removing temporary file...'
subprocess.call(['rm',pwd+'tmp.nc'])
