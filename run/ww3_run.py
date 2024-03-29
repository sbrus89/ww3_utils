import subprocess
import calendar
import pprint
import submission_script
import os
import glob
import sys
sys.dont_write_bytcode = True
import yaml
import pprint
import argparse
import time
import datetime

#####################################################################################################
#####################################################################################################

def get_date_ranges(date_start,date_end,days_per_run):

  date_range = []
  restart_interval = []
  date_run_start = datetime.datetime.strptime(date_start,'%m/%d/%Y')
  date_run_end = datetime.datetime.strptime(date_end,'%m/%d/%Y')
  run_length = datetime.timedelta(days=days_per_run) 
  keep_going = True

  date_range_start = date_run_start
  while keep_going:
 
    # Keep track date ranges 
    date_range_end = date_range_start + run_length

    # Handle the last date range 
    if date_range_end >= date_run_end:
      date_range_end = date_run_end 
      keep_going = False
      if date_range_end == date_range_start:
        break

    # Format dates
    date1 = datetime.datetime.strftime(date_range_start, '%Y%m%d')
    date2 = datetime.datetime.strftime(date_range_end,   '%Y%m%d')

    # Add date range to list
    date_range.append([date1 + ' 000000', 
                       date2 + ' 000000'])
    delta = date_range_end-date_range_start
    restart_interval.append(delta.days)
  
   # Update for next iteration
    date_range_start = date_range_end
  
  
  pprint.pprint(date_range)
  pprint.pprint(restart_interval)

  run_range = date_run_end - date_run_start
  if sum(restart_interval) != run_range.days:
    print("Error in date ranges")

  return date_range,restart_interval

#####################################################################################################
#####################################################################################################

def prep_grid(pwd,skip_existing=False,skip_grid=False,ignore_existing=False,require_existing=False):
 
  # Decide whether to prep grid
  if skip_grid:
     return
  if os.path.exists(pwd+'/mod_def.ww3'):
    if skip_existing:
      run_grid = 'n'
      print('mod_def.ww3 file exists, skipping ww3_grid')
    elif ignore_existing:
      run_grid = 'y'
    else:
      run_grid = input('mod_def.ww3 file exists, run ww3_grid? ')
  else:
    if require_existing:
      print('mod_def.ww3 file does not exist, quitting...')
      raise SystemExit(0)
    else:
      run_grid = 'y'

  # Create ww3_grid.inp and run ww3_grid
  if run_grid == 'y':
    subprocess.call(['python','ww3_grid.py'])
    subprocess.call([pwd+'/ww3_grid'])

#####################################################################################################
#####################################################################################################

def prep_ic(pwd,skip_existing=False,skip_strt=False,ignore_existing=False,require_existing=False):

  # Decide wheter to prep ic
  if skip_strt:
    return
  if os.path.lexists(pwd+'/restart.ww3'):
    if skip_existing:
      run_strt = 'n'
      print('restart.ww3 file exists, skipping ww3_strt')
    elif ignore_existing:
      run_strt = 'y'
    else:
      run_strt = input('restart.ww3 file exists, run ww3_strt? ')
  else:
    if require_existing:
      print('restart.ww3 file does not exist, quitting...')
      raise SystemExit(0)
    else:
      run_strt = 'y'

  # Create ww3_strt.inp and run ww3_strt
  if run_strt == 'y':
    subprocess.call(['python','ww3_strt.py'])
    subprocess.call(['srun','-n','4',pwd+'/ww3_strt'])

#####################################################################################################
#####################################################################################################

def prep_forcing(pwd,cfg,skip_existing=False,skip_prnc=False,ignore_existing=False,require_existing=False):

  if skip_prnc:
    return

  # Set forcings to default vaules if not present in config file
  if 'wind' not in cfg:
    cfg['wind'] = 'T'        # Always use wind forcing
  else:
    cfg['wind'] = 'T'
  if 'currents' not in cfg:  # Turn others off if not specified
    cfg['currents'] = 'F'
  if 'ssh' not in cfg:
    cfg['ssh'] = 'F'
  if 'ice' not in cfg:
    cfg['ice'] = 'F'
  if 'icebergs' not in cfg:
    cfg['icebergs'] = 'F'

  # Determine which forcings are requested
  forcings = []
  if cfg['wind'] == 'T':
    forcings.append('wind')
  if cfg['currents'] == 'T':
    forcings.append('currents')
  if cfg['ssh'] == 'T':
    forcings.append('ssh')
  if cfg['ice'] == 'T':
    forcings.append('ice')
  if cfg['icebergs'] == 'T':
    forcings.append('icebergs')
   
  ww3_files = {'wind':'wind.ww3','currents':'current.ww3','ssh':'level.ww3','ice':'ice.ww3','icebergs':'ice.ww3'}

  for forcing in forcings:
    if cfg[forcing] == 'T':

      # Decide whether to prep forcing
      if os.path.exists(pwd+'/'+ww3_files[forcing]):
        if skip_existing:
          run_prnc = 'n'
          print(forcing+'.ww3 file exisits, skipping ww3_prnc')
        elif ignore_existing:
          run_prnc = 'y'
        else:
          run_prnc = input(forcing+'.ww3 file exists, run ww3_prnc? ')
      else:
        if require_existing:
          print(forcing+'.ww3 file does not exist, quitting...')
          raise SystemExit(0)
        else:
          run_prnc = 'y'

      # Create ww3_prnc.inp, link correct data file, and run ww3_prnc
      if run_prnc == 'y':
 
        # if only one forcing is requested don't require the config file to be named with a '.forcing' extension
        # otherwise, copy the file with the '.forcing' file without the extension so ww3_prnc.py can read it (delete it later)
        if len(forcings) == 1 and  os.path.exists('ww3_prnc.config'):
          rm_config = False 
        else:
          subprocess.call(['cp','ww3_prnc.config.'+forcing,'ww3_prnc.config'])
          rm_config = True

        subprocess.call(['python','ww3_prnc.py'])
        direc = cfg[forcing+"_direc"]
        ncfile = glob.glob(direc+'*ww3.nc')[0]
        subprocess.call(['ln','-sf',ncfile,forcing+'.nc'])
        subprocess.call(['srun','-n','36',pwd+'/ww3_prnc'])
        if rm_config:
          subprocess.call(['rm','ww3_prnc.config'])     

#####################################################################################################
#####################################################################################################

def prep_shel(pwd,skip_existing=False,ignore_existing=False,require_existing=False):

  # Decide whether to prep shel
  if os.path.exists(pwd+'/ww3_shel.inp'):
    if skip_existing:
      gen_shel = 'n'
      print('ww3_shel.inp file exists, skipping ww3_shel.py')
    elif ignore_existing:
      gen_shel = 'y'
    else:
      gen_shel = input('ww3_shel.inp file exists, run ww3_shel.py?')
  else:
    if require_existing:
      print('ww3_shel.inp does not exist, quitting...')
      raise SystemExit(0)
    else:
      gen_shel = 'y'

  # Create ww3_shel.inp file, first updating the start/end times in ww3_shel.config
  if gen_shel == 'y':
    lines = open(pwd+'/ww3_shel.config').read().splitlines()
    for i,line in enumerate(lines):
      if line and line.split()[0] == 'start_time':
         start_time = datetime.datetime.strptime(cfg["date_start"],'%m/%d/%Y')
         lines[i] = "start_time     : '"+datetime.datetime.strftime(start_time,'%Y%m%d')+" 000000'"
      if line and line.split()[0] == 'end_time':
         end_time = datetime.datetime.strptime(cfg["date_end"],'%m/%d/%Y')
         lines[i] = "end_time       : '"+datetime.datetime.strftime(end_time,'%Y%m%d')+" 000000'"
    f = open(pwd+'/ww3_shel.config','w')
    f.write('\n'.join(lines))
    f.close()
    subprocess.call(['python','ww3_shel.py'])

#####################################################################################################
#####################################################################################################

def restart_commands(i,restart_interval,date_range):

  start = date_range[i][0]
  end = date_range[i][1]

  # Setup initial restart  
  if i == 0:
    interval = str(int(restart_interval[i])*24*3600)
    pre_cmds = ['python ww3_restart.py --restart_time="'+start+'" '+ \
                                       '--stop_time="'+end+'" '+ \
                                       '--restart_interval="'+interval+'" '+ \
                                       '--skip_log']
  else:
    pre_cmds = []

  # Setup next restart
  if i+1 < len(date_range):
    restart = date_range[i+1][0]
    restart_end = date_range[i+1][1]
    interval = str(int(restart_interval[i+1])*24*3600)

    post_cmds = ['python ww3_restart.py --restart_time="'+restart+'" '+ \
                                       '--stop_time="'+restart_end+'" '+ \
                                       '--restart_interval="'+interval+'"']
  # Rename last output files
  else:
    post_cmds = ['python ww3_restart.py']

  return pre_cmds,post_cmds

#####################################################################################################
#####################################################################################################

if __name__ == '__main__':

  parser = argparse.ArgumentParser()
  parser.add_argument('--submit'          , action='store_true', help='submit the jobs with dependencies')
  parser.add_argument('--test'            , action='store_true', help='write out test files to test restart setup')
  parser.add_argument('--skip_existing'   , action='store_true', help='skip over actions for existing file')
  parser.add_argument('--ignore_existing' , action='store_true', help='re-do actions for existing file')
  parser.add_argument('--require_existing', action='store_true', help='require files to elready xist')
  parser.add_argument('--skip_grid'       , action='store_true', help='skip running ww3_grid')
  parser.add_argument('--skip_strt'       , action='store_true', help='skip running ww3_strt')
  parser.add_argument('--skip_prnc'       , action='store_true', help='skip running ww3_prnc')
  args = parser.parse_args()

  pwd = os.getcwd()
 
  # Read configuration file
  f = open(pwd+'/ww3_run.config')
  cfg = yaml.load(f, Loader=yaml.Loader)
  pprint.pprint(cfg)


  # Create ww3_grid.inp and run ww3_grid
  prep_grid(pwd, args.skip_existing, args.skip_grid, args.ignore_existing)
  
  # Create ww3_strt.inp and run ww3_strt
  prep_ic(pwd, args.skip_existing, args.skip_strt, args.ignore_existing)

  # Create ww3_prnc.inp, link correct data file, and run ww3_prnc
  prep_forcing(pwd, cfg, args.skip_existing, args.skip_prnc, args.ignore_existing)

  # Create ww3_shel.inp file, first updating the start/end times in ww3_shel.config
  prep_shel(pwd, args.skip_existing, args.ignore_existing)


  # Get the start and end dates for each run
  date_range,restart_interval = get_date_ranges(cfg["date_start"],cfg["date_end"],cfg["days_per_run"])

  # Check for existing restart files
  existing_restarts = glob.glob(pwd+'/restart.ww3.*')  
  

  # Setup the submission scripts for the series of runs
  sub_count = 0
  for i in range(len(date_range)):
    start = date_range[i][0]
    end = date_range[i][1]

    sub_file = pwd+'/ww3.'+start.replace(' ','_')+'-'+end.replace(' ','_')+'.sub'

 
    # Set the estimated runtime
    if 'runtime' in cfg:
      runtime = cfg['runtime']
    elif 'runtime_per_day' in cfg:
      runtime = restart_interval[i]*cfg['runtime_per_day'] 
    else:
      runtime = None

    pre_cmds,post_cmds = restart_commands(i,restart_interval,date_range) 

    # Write the submission script
    submission_script.write_submission_script(machine=cfg["machine"],
                                              ncores=cfg["ncores"],
                                              job_name=cfg["job_name"],
                                              queue=cfg["queue"],
                                              exe=cfg["exe"],
                                              time=runtime,
                                              filename=sub_file,
                                              pre_cmds=pre_cmds,
                                              post_cmds=post_cmds)

    # Check if end restart files exist
    restart_end   = pwd+'/restart.ww3.'+end.replace(' ','_')
    if restart_end in existing_restarts:
      skip_submit = True
    else:
      skip_submit = False

    # Submit the runs with depenencies
    if args.submit:
      if os.path.isfile(pwd+'/'+cfg['exe']) and not skip_submit :
        if sub_count > 0:
          run_cmd = ['sbatch','--dependency=afterok:'+job_id,sub_file]
        else:
          run_cmd = ['sbatch',sub_file]
  
        sub_count = sub_count + 1

        print(' '.join(run_cmd))
        output = subprocess.Popen(run_cmd, stdout=subprocess.PIPE).communicate()[0]
        print(output)
  
        job_id = output.split()[-1]
        
        time.sleep(3)
    
    #-------------------------------------------------------------------------------
    # Run test on dummy files
    if args.test:
      if len(pre_cmds) > 0:
        subprocess.call(pre_cmds[0],shell=True)

      pause = input("Beginning model run...")

      f = open('restart001.ww3','w')
      f.write('Test restart file for '+end)
      f.close()

      f = open('out_grd.ww3','w')
      f.write('Test field output file for '+start+'-'+end)
      f.close()

      f = open('out_pnt.ww3','w')
      f.write('Test point output file for '+start+'-'+end)
      f.close()

      f = open('log.ww3','w')
      f.write('Test log file for '+start+'-'+end+'\n')
      f.write('                                        |       input       |     output    |'+'\n')
      f.write('      step | pass |    date      time   | b w l c i i1 i5 d | g p t r b f c |'+'\n')
      f.write('   --------|------|---------------------|-------------------|---------------|'+'\n')

      date = start.split()[0]
      time = start.split()[1]
      start_date = date[0:4]+'/'+date[4:6]+'/'+date[6:]
      start_time = time[0:2]+':'+time[2:4]+':'+time[4:]
      f.write('       0   |   0  | '+start_date+' '+start_time+' |   X               | X X           |'+'\n')

      date = end.split()[0]
      time = end.split()[1]
      end_date = date[0:4]+'/'+date[4:6]+'/'+date[6:]
      end_time = time[0:2]+':'+time[2:4]+':'+time[4:]
      f.write('       0   |   0  | '+end_date+' '+end_time+    ' |   X               | X X   X       |'+'\n')
      f.write('   --------+------+---------------------+-------------------+---------------+'+'\n')
      f.close()

      pause = input("End of model run")
      subprocess.call(post_cmds[0],shell=True)
      pause = input("End of restart setup")
     
    
 
    
