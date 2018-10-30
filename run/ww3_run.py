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




def get_date_ranges(year,month,days_per_run):

  days_in_month = calendar.monthrange(year,month)[1]
  
  date_range = []
  restart_interval = []
  date_range_start = 1
  keep_going = True

  while keep_going:
 
    # Keep track date ranges 
    date_range_end = date_range_start + days_per_run
  
    # Format dates 
    date1 = str(date_range_start).zfill(2)
    date2 = str(date_range_end).zfill(2)
    month1 = str(month).zfill(2)
    month2 = month1

    # Handle the last date range (ends at beginning of the first day of the next month) 
    if date_range_end > days_in_month:
      month2 = str(month + 1).zfill(2)
      date2 = '01'  
      date_range_end = days_in_month+1
      keep_going = False

    # Add date range to list
    date_range.append([str(year)+month1+date1 + ' 000000', 
                       str(year)+month2+date2 + ' 000000'])
    restart_interval.append(date_range_end-date_range_start)
  
   # Update for next iteration
    date_range_start = date_range_end
  
  
  print days_in_month
  pprint.pprint(date_range)
  pprint.pprint(restart_interval)

  if sum(restart_interval) != days_in_month:
    print "Error in date ranges"

  return date_range,restart_interval

#####################################################################################################
#####################################################################################################

if __name__ == '__main__':

  parser = argparse.ArgumentParser()
  parser.add_argument('--submit', action='store_true', help='submit the jobs with dependencies')
  parser.add_argument('--test'  , action='store_true', help='write out test files to test restart setup')
  args = parser.parse_args()

  pwd = os.getcwd()
 
  # Read configuration file
  f = open(pwd+'/ww3_run.config')
  cfg = yaml.load(f)
  pprint.pprint(cfg)


  # Create ww3_grid.inp and run ww3_grid
  run_grid = 'y'
  if os.path.exists(pwd+'/mod_def.ww3'):
    run_grid = raw_input('mod_def.ww3 file exists, run ww3_grid? ')
  if run_grid == 'y':
    subprocess.call(['python','ww3_grid.py'])
    subprocess.call([pwd+'/ww3_grid'])
  
  # Create ww3_strt.inp and run ww3_strt
  run_strt = 'y'
  if os.path.exists(pwd+'/restart.ww3'):
    run_strt = raw_input('restart.ww3 file exists, run ww3_strt? ')
  if run_strt == 'y':
    subprocess.call(['python','ww3_strt.py'])
    subprocess.call([pwd+'/ww3_strt'])

  # Set forcings to default vaules if not present in config file
  if 'wind' not in cfg:
    cfg['wind'] = 'T'        # Always use wind forcing
  else:
    cfg['wind'] = 'T'
  if 'currents' not in cfg:  # Turn others off if not specified
    cfg['currents'] = 'F'
  if 'ssh' not in cfg:
    cfg['ssh'] = 'F'

  # Determine which forcings are requested
  forcings = []
  if cfg['wind'] == 'T':
    forcings.append('wind')
  if cfg['currents'] == 'T':
    forcings.append('currents')
  if cfg['ssh'] == 'T':
    forcings.append('ssh')
   
  ww3_files = {'wind':'wind.ww3','currents':'current.ww3','ssh':'level.ww3'}

  # Create ww3_prnc.inp, link correct data file, and run ww3_prnc
  for forcing in forcings:
    if cfg[forcing] == 'T':
      run_prnc = 'y'
      if os.path.exists(pwd+'/'+ww3_files[forcing]):
        run_prnc = raw_input(forcing+'.ww3 file exists, run ww3_prnc? ')
      if run_prnc == 'y':
 
        # if only one forcing is requested don't require the config file to be named with a '.forcing' extension
        # otherwise, copy the file with the '.forcing' file without the extension so ww3_prnc.py can read it (delete it later)
        if len(forcings) == 1 and  os.path.exists('ww3_prnc.config'):
          rm_config = False 
        else:
          subprocess.call(['cp','ww3_prnc.config.'+forcing,'ww3_prnc.config'])
          rm_config = True

        subprocess.call(['python','ww3_prnc.py'])
        month_name = calendar.month_name[cfg["month"]].lower()
        direc = cfg[forcing+"_direc"]+month_name+'/'
        ncfile = glob.glob(direc+'*ww3.nc')[0]
        subprocess.call(['ln','-sf',ncfile,forcing+'.nc'])
        subprocess.call([pwd+'/ww3_prnc'])
        if rm_config:
          subprocess.call(['rm','ww3_prnc.config'])     

  # Create ww3_shel.inp file, first updating the start/end times in ww3_shel.config
  gen_shel = 'y'
  if os.path.exists(pwd+'/ww3_shel.inp'):
    gen_shel = raw_input('ww3_shel.inp file exists, run ww3_shel.py?')
  if gen_shel == 'y':
    lines = open(pwd+'/ww3_shel.config').read().splitlines()
    for i,line in enumerate(lines):
      if line and line.split()[0] == 'start_time':
         lines[i] = "start_time     : '"+str(cfg["year"])+str(cfg["month"]).zfill(2)+"01 000000'"
      if line and line.split()[0] == 'end_time':
         lines[i] = "end_time       : '"+str(cfg["year"])+str(cfg["month"]+1).zfill(2)+"01 000000'"
    f = open(pwd+'/ww3_shel.config','w')
    f.write('\n'.join(lines))
    f.close()
    subprocess.call(['python','ww3_shel.py'])


  date_range,restart_interval = get_date_ranges(cfg["year"],cfg["month"],cfg["days_per_run"])

  for i in range(len(date_range)):
    start = date_range[i][0]
    end = date_range[i][1]

    sub_file = pwd+'/ww3.'+start.replace(' ','_')+'-'+end.replace(' ','_')+'.sub'

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
 

    # Write the submission script
    submission_script.write_submission_script(machine=cfg["machine"],
                                              ncores=cfg["ncores"],
                                              job_name=cfg["job_name"],
                                              queue=cfg["queue"],
                                              exe=cfg["exe"],
                                              filename=sub_file,
                                              pre_cmds=pre_cmds,
                                              post_cmds=post_cmds)


    # Submit the runs with depenencies
    if args.submit:
      if os.path.isfile(pwd+'/'+cfg['exe']):
        if i > 0:
          run_cmd = ['sbatch','--dependency=afterok:'+job_id,sub_file]
        else:
          run_cmd = ['sbatch',sub_file]
  
        print ' '.join(run_cmd)
        output = subprocess.Popen(run_cmd, stdout=subprocess.PIPE).communicate()[0]
        print output
  
        job_id = output.split()[-1]
        
        time.sleep(3)
    
    #-------------------------------------------------------------------------------
    # Run test on dummy files
    if args.test:
      if len(pre_cmds) > 0:
        subprocess.call(pre_cmds[0],shell=True)

      pause = raw_input("Beginning model run...")

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

      pause = raw_input("End of model run")
      subprocess.call(post_cmds[0],shell=True)
      pause = raw_input("End of restart setup")
     
    
 
    
