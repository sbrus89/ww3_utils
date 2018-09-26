import subprocess
import os

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

pwd = os.getcwd()+'/'

first_file  = 'wnd10mx0.5.gdas.200505.grb2.nc' 
second_file = 'wnd10mx0.5.gdas.200506.grb2.nc'

print 'Concatenatng files...'
subprocess.call(['ncrcat',pwd+first_file,second_file,pwd+'tmp.nc'])
print 'Fixing latitude dimension...'
subprocess.call(['ncpdq','-O','-a','-lat',pwd+'tmp.nc',pwd+second_file.replace('grb2','ww3')])
print 'Removing temporary file...'
subprocess.call(['rm',pwd+'tmp.nc'])
