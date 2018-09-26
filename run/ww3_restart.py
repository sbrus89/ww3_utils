import os
import glob
import subprocess

#############################################################################################################

def link_restart(restart_time,start_time):
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
  last_restart = 'restart'+str(restart_max).zfill(3)+'.ww3'

  # If restartXXX.ww3 files were found, link restart.ww3 to the most recent one
  if restart_max > 0:

    # Save initial restart.ww3
    if restart_exists == True and restart_link == False:                    
      subprocess.call(['mv', 'restart.ww3', 'restart.ww3.'+start_time])     

    # Rename most recent restart file and link to restart.ww3
    subprocess.call(['mv',last_restart,'restart.ww3.'+restart_time])
    subprocess.call(['ln','-sf','restart.ww3.'+restart_time,'restart.ww3'])

#############################################################################################################

def restart_time():
  pwd = os.getcwd()

  log_file = pwd+'/log.ww3'

  if not os.path.isfile(log_file):
    print 'log.ww3 file does not exist'
    raise SystemExit(0)

  f = open(log_file,'r')
  lines = f.read().splitlines()
 
  output_intervals = False
  restart_time = ''
  start_time = ''
  for line in lines:

    # Find lines inside output table
    if line.find('--------|------|---------------------|-------------------|---------------|') >= 0:
      output_intervals = True
      continue
    if line.find('--------+------+---------------------+-------------------+---------------+') >= 0:
      output_intervals = False

    # Process lines inside output table
    if output_intervals == True:
      cols = line.split('|')
      
      # Find date and time
      date_time = cols[2].strip().split()
      if len(date_time) > 1:
        date = date_time[0]
        time = date_time[1]
      else:
        time = date_time[0]
      date_time_formatted = ''.join(date.split('/')) + ' ' + ''.join(time.split(':'))

      # Find start time
      if start_time == '':
        start_time = date_time_formatted

      # Find restart output
      output = cols[4]
      if output[7] == 'X' or output[7] == 'L':
        restart_time = date_time_formatted

  return restart_time,start_time

#############################################################################################################

def rename_outputs(restart_time,start_time):

  pwd = os.getcwd()

  files = ['out_grd.ww3','out_pnt.ww3']
  interval = start_time + '-' + restart_time

  for name in files:
    f = pwd + '/' + name
    if os.path.isfile(f):
      subprocess.call(['mv',f,f+'.'+interval])

#############################################################################################################

def update_shel_input(restart_time):

  pwd = os.getcwd()

  shel_inp = pwd + '/ww3_shel.inp'

  f = open(shel_inp,'r')
  lines = f.read().splitlines()
  f.close()

  i = 0
  for j,line in enumerate(lines):
    if line[0] != '$':
      i = i + 1
      if i == 8:
        start_time_line = j
        break

  lines[start_time_line] = '   '+restart_time

  f = open(shel_inp,'w')
  for line in lines:
    f.write(line + '\n')
  f.close()

#############################################################################################################
 
if __name__ == '__main__':
  restart_time,start_time = restart_time()
  update_shel_input(restart_time)
  restart_time = restart_time.replace(' ','_')
  start_time = start_time.replace(' ','_')
  link_restart(restart_time,start_time)
  rename_outputs(restart_time,start_time)
