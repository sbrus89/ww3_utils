import os
import glob
import subprocess

def setup_restart():
  pwd = os.getcwd()
  
  restart_file = pwd+'/restart.ww3'
  
  # Check if restart.ww3 exists
  restart_exists = False
  if os.path.isfile(restart_file):
    restart_exists = True
  
  # Check if restart.ww3 is a link
  restart_link = False
  if os.path.islink(restart_file):
    restart_link = True
  
  # Look for most recent restartXXX.ww3 file
  restart_max = 0
  restart_files = glob.glob('restart*.ww3')
  for f in restart_files:
    print f
    restart_id = f[ f.find('restart')+len('restart'):f.rfind('.')]
    if restart_id.strip():
      restart_n = int(restart_id)
      if restart_n > restart_max:
        restart_max = restart_n
  
  # If restartXXX.ww3 files were found, link restart.ww3 to the most recent one
  if restart_max > 0:
    if restart_exists == True and restart_link == False:
      subprocess.call(['mv', 'restart.ww3', 'restart000.ww3'])
    subprocess.call(['ln','-sf','restart'+str(restart_max).zfill(3)+'.ww3','restart.ww3'])
 
if __name__ == '__main__':
  setup_restart() 
