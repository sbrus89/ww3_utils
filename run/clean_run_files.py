import os
import yaml
import glob
import pprint 
import calendar
import subprocess
import datetime

pwd = os.getcwd()

# Read ww3_run.config to get month
f = open(pwd+'/ww3_run.config')
cfg = yaml.load(f)
pprint.pprint(cfg)

# Define some date strings for later
start_date = datetime.datetime.strptime(cfg["date_start"],'%m/%d/%Y')
start_date = datetime.datetime.strftime(start_date,'%Y%m%d')
end_date   = datetime.datetime.strptime(cfg["date_end"],  '%m/%d/%Y')
end_date   = datetime.datetime.strftime(end_date,'%Y%m%d')
year = start_date[0:4]

# Get names of restart files
first_restart = 'restart.ww3.'+start_date+'_000000'
last_restart  = 'restart.ww3.'+end_date  +'_000000'
all_restarts = glob.glob('restart.ww3.*') 

# Get list of intermediate restart files
if first_restart not in all_restarts:
  print 'first restart file not found'
  keep_going = raw_input('continue? ')
  if keep_going != 'y':
    raise SystemExit(0)
else:
  all_restarts.remove(first_restart)
if last_restart not in all_restarts:
  print 'last restart file not found'
  keep_going = raw_input('continue? ')
  if keep_going != 'y':
    raise SystemExit(0)
else:
  all_restarts.remove(last_restart)
other_restarts = all_restarts

# List of commands to be run
tasks = [{'cmd':'mv', 'files':'log*',       'dest':'results/model_output/log_files'},
         {'cmd':'mv', 'files':'out_pnt*',   'dest':'results/model_output/points'},
         {'cmd':'mv', 'files':'out_grd*',   'dest':'results/model_output/fields'},
         {'cmd':'mv', 'files':'ww3_*.o*',   'dest':'results/model_output/screen_output'},
         {'cmd':'mv', 'files':'ww3_*.e*',   'dest':'results/model_output/screen_output'},
         {'cmd':'mv', 'files':first_restart,'dest':'results/model_output/restarts'},
         {'cmd':'mv', 'files':'*.sub',      'dest':'results/submission_scripts'},
         {'cmd':'rm', 'files':other_restarts}]


for task in tasks:

  # Create destination directory, if necessary
  if 'dest' in task:
    dest_direc = pwd+'/'+year+'/'+task['dest']
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
    command.extend(files)
    if 'dest' in task:
      command.append(dest_direc)

    # Run command
    subprocess.call(command)
