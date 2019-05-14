import subprocess
import os
import glob

# Combines CFSR NetCDF files for a given time period
# into a single NetCDF file that can be read by ww3_prnc.
#
# For wind data, the first timesnap for the time period is
# contained in its own file:
#
# The first file contains only the 00 hour data for the first 
#    day of the month contained in the second file.
# The second file contains the rest of the hourly data for that month, 
#    beginning with the 01 hour of the first day.
#
# The latitude dimension also has to be altered for ww3_prnc
#
# For ocean current data, the u and v components must first be
# appended into a single file.



pwd = os.getcwd()+'/'

filetypes = ['icecon','prmsl','wnd10m','ocnu']


for ft in filetypes:

  files = sorted(glob.glob(pwd+ft+'*.nc'))
  
  if len(files) < 1:
    continue

  filenames = [f.split('/')[-1] for f in files]
  print filenames

  date1 = files[1].split('.')[-3]
  date2 = files[-1].split('.')[-3]
  
  if ft == 'ocnu':
    print 'Appending files...'  
    for f in filenames:
      fileA = f
      fileB = f.replace('ocnu','ocnv')
      subprocess.call(['ncks','-A',pwd+fileB,pwd+fileA])
  
  print 'Concatenatng files...'
  command = ['ncrcat']
  command.extend(filenames)
  command.append(pwd+'tmp.nc')
  subprocess.call(command)
  print 'Fixing latitude dimension...'
  subprocess.call(['ncpdq','-O','-a','-lat',pwd+'tmp.nc',pwd+ft+'.'+date1+'-'+date2+'.ww3.nc'])
  print 'Removing temporary file...'
  subprocess.call(['rm',pwd+'tmp.nc'])
