import subprocess
import calendar
import pprint
import submission_script
import os

year = 2005
month = 9 
days_per_run = 5

test = False 

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

  pwd = os.getcwd()

  date_range,restart_interval = get_date_ranges(year,month,days_per_run)

  for i in range(len(date_range)):
    start = date_range[i][0]
    end = date_range[i][1]

    sub_file = pwd+'/ww3.'+start.replace(' ','_')+'-'+end.replace(' ','_')+'.sub'

    # Setup initial restart  
    if i == 0:
      interval = str(int(restart_interval[i])*24*3600)
      pre_cmds = ['python ww3_restart.py --restart_time="'+start+'" '+ \
                                         '--stop_time="'+end+'" '+ \
                                         '--restart_interval="'+interval+'"']
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
    submission_script.write_submission_script(machine='grizzly',
                                              ncores=216,
                                              job_name='ww3_glo_15m',
                                              queue='interactive',
                                              exe='ww3_shel',
                                              filename=sub_file,
                                              pre_cmds=pre_cmds,
                                              post_cmds=post_cmds)


    # Submit the runs with depenencies
    if os.path.isfile(pwd+'/ww3_shel'):
      if i > 0:
        run_cmd = ['sbatch','--dependency=afterok:'+job_id,sub_file]
      else:
        run_cmd = ['sbatch',sub_file]

      print ' '.join(run_cmd)
      output = subprocess.Popen(run_cmd, stdout=subprocess.PIPE).communicate()[0]
      print output

      job_id = output.split()[-1]

    
    #-------------------------------------------------------------------------------
    # Run test on dummy files
    if test:
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
     
    
 
    
