import os
import yaml
import glob
import pprint 
import calendar
import subprocess

pwd = os.getcwd()

# Read ww3_run.config to get month
f = open(pwd+'/ww3_run.config')
cfg = yaml.load(f)
pprint.pprint(cfg)

# Define some strings for later
restart_file = 'restart.ww3.'+str(cfg['year'])+str(cfg['month']).zfill(2)+'01_000000'
other_restarts = restart_file.strip('01_000000')+'*'
month_name = calendar.month_name[cfg["month"]].lower()

# List of commands to be run
tasks = [{'cmd':'mv', 'files':'log*',      'dest':'model_output/log_files'},
         {'cmd':'mv', 'files':'out_pnt*',  'dest':'model_output/points'},
         {'cmd':'mv', 'files':'out_grd*',  'dest':'model_output/fields'},
         {'cmd':'mv', 'files':'ww3_*.o*',  'dest':'model_output/screen_output'},
         {'cmd':'mv', 'files':'ww3_*.e*',  'dest':'model_output/screen_output'},
         {'cmd':'mv', 'files':restart_file,'dest':'model_output'},
         {'cmd':'mv', 'files':'*.sub',     'dest':'subission_scripts'},
         {'cmd':'rm', 'files':other_restarts}]


for task in tasks:

  # Create destination directory, if necessary
  if 'dest' in task:
    dest_direc = pwd+'/post_processing/'+month_name+'/'+task['dest']
    if not os.path.exists(dest_direc):
      subprocess.call(['mkdir','-p',dest_direc])

  # Find list of files
  files = glob.glob(task['files'])


  if len(files) > 0:

    # Build command
    command = [task['cmd']]
    command.extend(files)
    if 'dest' in task:
      command.append(dest_direc)

    # Run command
    subprocess.call(command)
