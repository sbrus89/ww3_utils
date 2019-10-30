import os
import yaml
import glob
import pprint 
import calendar
import subprocess
import datetime

pwd = os.getcwd()

# Get names of restart files
all_restarts = glob.glob('restart.ww3.*') 
if len(all_restarts) > 0:
  frmt = '%Y%m%d_%H%M%S'
  restart_dates = []
  for restart_file in all_restarts:
    date = restart_file.split('.')[-1]
    restart_dates.append(datetime.datetime.strptime(date,frmt))
  restart_dates.sort()
  start_date = datetime.datetime.strftime(restart_dates[0],frmt)
  end_date = datetime.datetime.strftime(restart_dates[-1],frmt)
  first_restart = 'restart.ww3.'+start_date
  last_restart  = 'restart.ww3.'+end_date  
  # Get list of intermediate restart files
  #all_restarts.remove(first_restart)
  all_restarts.remove(last_restart)
  other_restarts = all_restarts
else:
  first_restart = 'restart.ww3.XXXXXXXX_XXXXXX'
  last_restart  = 'restart.ww3.XXXXXXXX_XXXXXX'  
  other_restarts = [] 


# List of commands to be run
tasks = [{'cmd':'mv', 'opt':'' , 'files':'log*',       'dest':'results/model_output/log_files'},
         {'cmd':'mv', 'opt':'' ,'files':'out_pnt*',   'dest':'results/model_output/points'},
         {'cmd':'mv', 'opt':'' ,'files':'out_grd*',   'dest':'results/model_output/fields'},
         {'cmd':'mv', 'opt':'' ,'files':'run*.o*',    'dest':'results/model_output/screen_output'},
         {'cmd':'mv', 'opt':'' ,'files':'run*.e*',    'dest':'results/model_output/screen_output'},
         {'cmd':'mv', 'opt':'' ,'files':last_restart, 'dest':'results/model_output/restarts'},
         {'cmd':'mv', 'opt':'' ,'files':'*.sub',      'dest':'results/submission_scripts'},
         {'cmd':'rm', 'opt':'' ,'files':other_restarts},
         {'cmd':'ln', 'opt':'-s' ,'files':'results/model_output/restarts/restart.ww3.*'}]


for task in tasks:

  # Create destination directory, if necessary
  if 'dest' in task:
    dest_direc = pwd+'/'+task['dest']
    if not os.path.exists(dest_direc):
      subprocess.call(['mkdir','-p',dest_direc])

  # Find list of files
  try:
    files = glob.glob(task['files'])
  except:
    files = task['files']

  if len(files) > 0:

    # Build command
    command = [task['cmd']]
    command.append(task['opt'])
    command.extend(files)
    if 'dest' in task:
      command.append(dest_direc)

    # Run command
    print(command)
    subprocess.call(command)
