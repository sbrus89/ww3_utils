import subprocess
import os
import glob

pwd = os.getcwd()+'/'
 
filetypes = ['icecon','prmsl','wnd10m']

for ft in filetypes:

  files = sorted(glob.glob(pwd+ft+'*.nc'))
  print(files) 

  filenames = [f.split('/')[-1] for f in files]
  print(filenames)

  if len(filenames) > 0:
    command = ['ncrcat']
    command.extend(filenames)
    command.append(ft+'.nc')
    print(command)
    subprocess.call(command)
