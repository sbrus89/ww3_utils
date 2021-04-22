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

def find_restart_time(restart_time_expected=None):
  pwd = os.getcwd()

  log_file = pwd+'/log.ww3'

  if not os.path.isfile(log_file):
    print('log.ww3 file does not exist')
    raise SystemExit(0)

  f = open(log_file,'r')
  lines = f.read().splitlines()
 
  output_intervals = False
  restart_found = False
  restart_time = ''
  start_time = ''
  end_time = ''
  n = 0
  for i,line in enumerate(lines):

    n = n + 1

    # Find lines inside output table
    if line.find('--------|------|---------------------|-------------------|---------------|') >= 0:
      output_intervals = True
      n = 0 
      continue
    if line.find('--------+------+---------------------+-------------------+---------------+') >= 0 and n > 1: # For whatever reason, for E3SM output this line is  
      output_intervals = False                                                                                 # right afer the ----|----|----| line with output 
    elif n == 1:                                                                                               # following it. This prevents if from being flaged
      continue                                                                                                 # as the end of the file    

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
        restart_found = True

  # Set restart_time as the final time if no restart was written
  if restart_found == False:
    end_time = date_time_formatted 
    print("no restart file written")
  else:
    if restart_time_expected and restart_time != restart_time_expected:
      print("restart time is different from what is expected")
      raise SystemExit(0)

  return restart_time,start_time,end_time

#############################################################################################################

def rename_outputs(restart_time,start_time):

  pwd = os.getcwd()

  files = ['out_grd.ww3','out_pnt.ww3','log.ww3']
  interval = start_time + '-' + restart_time

  for name in files:
    f = pwd + '/' + name
    if os.path.isfile(f):
      subprocess.call(['mv',f,f+'.'+interval])

#############################################################################################################

def update_shel_input(restart_time,stop_time=None,restart_interval=None):

  pwd = os.getcwd()

  shel_inp = pwd + '/ww3_shel.inp'
  if not os.path.exists(shel_inp):
    return

  print("Updating: "+shel_inp)

  f = open(shel_inp,'r')
  lines = f.read().splitlines()
  f.close()

  replace_ww3_shel_inp_line(lines,8,restart_time)
  if stop_time:
    replace_ww3_shel_inp_line(lines,9,stop_time)
  if restart_interval:
    replace_ww3_shel_inp_line(lines,17,'  '.join([restart_time,restart_interval,stop_time]))

  f = open(shel_inp,'w')
  for line in lines:
    f.write(line + '\n')
  f.close()

#############################################################################################################

def replace_ww3_shel_inp_line(file_lines,nline,new_text):

  # nline number represents requried, non-comment lines (skips all station locations)
  # 1-7: define input to be used
  # 8: start time
  # 9: stop time
  # 10: define output data
  # 11-13: mean wave parameters output
  # 14: point output
  # 15: 'STOPSTING'
  # 16: track output
  # 17: restart output
  # 18: boundary output
  # 19: separated output wave field

  stop_string = False
  i = 0
  for j,line in enumerate(file_lines):
    if line[0] != '$':
      if i > 14 and stop_string == False:
        if line.find('STOPSTRING') > 0:
          stop_string = True
      else:
        i = i + 1

      if i == nline:
        input_line = j
        break

  file_lines[input_line] = '   '+new_text


#############################################################################################################

def setup_restart(restart_time=None,stop_time=None,restart_interval=None,check_log=True):

  if check_log:
    # Get the most recent restart time from the log.ww3 file, If a restart time is 
    # provided, check to make sure it agrees with the log file.
    restart_time,start_time,end_time = find_restart_time(restart_time)

  if restart_time:
    update_shel_input(restart_time,stop_time,restart_interval)
    
    if check_log:
      restart_time = restart_time.replace(' ','_')
      start_time = start_time.replace(' ','_')

      link_restart(restart_time,start_time)

      rename_outputs(restart_time,start_time)
  else:
    end_time = end_time.replace(' ','_')
    start_time = start_time.replace(' ','_')

    rename_outputs(end_time,start_time)

#############################################################################################################

if __name__ == '__main__':
  import argparse

  parser = argparse.ArgumentParser()
  parser.add_argument("--restart_time"     , type=str,  help="time for model restart (YYYYMMDD hhmmss)")
  parser.add_argument("--stop_time"        , type=str,  help="time for model end     (YYYYMMDD hhmmss)")
  parser.add_argument("--restart_interval" , type=str,  help="interval to write restartfile (seconds)")   
  parser.add_argument("--skip_log"         , action='store_false', help="don't check the log.ww3 file for the latest restart time")
  args = parser.parse_args()

  setup_restart(args.restart_time, args.stop_time, args.restart_interval, args.skip_log)




