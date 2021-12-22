import subprocess
import os
import glob
import yaml
import pprint
import ww3_log
import datetime

pwd = os.getcwd()

# Output type options
out_types = {'wave':{'type':'2','prefix':'ww3.','subtype':'2'},  # mean wave parameters 
             'met': {'type':'2','prefix':'cfsr.','subtype':'1'}, # depth, current, wind
             'spec':{'type':'1','prefix':'spec.','subtype':'3'}} 
###################################################################################################
###################################################################################################

def replace_ww3_ounp_inp_line(comment,opt1,opt2=None,opt3=None,opt4=None):

  # Find the requested line (nline) in the ww3_ounp input file
  founp = open(pwd+'/ww3_ounp.inp','r')
  lines = founp.read().splitlines()
  for n,line in enumerate(lines):
    if line.find(comment) > 0:
      break
  n = n+1
  
  # Replace the line with the new information (opt1 and opt2) 
  line_info = lines[n].split()
  line_info[0] = opt1
  if opt2:
    line_info[1] = opt2
  if opt3:
    line_info[2] = opt3
  if opt4:
    line_info[3] = opt4
  lines[n] = '   '+'  '.join(line_info)
  founp.close()
  
  # Re-write the ww3_ounp input file
  founp = open(pwd+'/ww3_ounp.inp','w')
  founp.write('\n'.join(lines))
  founp.close()

###################################################################################################
###################################################################################################

if __name__ == '__main__':

  f = open(pwd+'/process_points.config')
  cfg = yaml.load(f,yaml.Loader)
  pprint.pprint(cfg)

  # Check if the ww3_ounp input file exists
  if not os.path.isfile(pwd+'/ww3_ounp.inp'):
    print('ww3_ounp.inp not found')
    raise SystemExit(0)
  
  # Link the mod_def.ww3 file to the current directory
  subprocess.call(['ln','-sf',cfg['run_direc']+'mod_def.ww3',pwd])
  
  # Loop over all out_pnt.ww3.YYYYMMDD_HHMMSS-YYMMDD_HHMMSS files
  pnt_files = sorted(glob.glob(cfg['output_direc']+'out_pnt.ww3*'))
  log_files = sorted(glob.glob(cfg['log_direc']+'log.ww3*'))
  for i in range(len(pnt_files)):
  
    f = pnt_files[i]
    # Link the out_pnt.ww3 file to the current directory
    subprocess.call(['ln','-sf',f,pwd+'/out_pnt.ww3'])
  
    # Find the start and end dates from the filename
    #date_range = f.split(".")[-1]  
    #start_date_time = date_range.split('-')[0]
    #start_date = start_date_time.split('_')[0]
    #start_time = start_date_time.split('_')[1]
    #print(start_date,start_time)

    # Find output interval and number of outputs
    restart_output_times,gridded_output_times,point_output_times,start,end = ww3_log.find_output_times(log_files[i])
    noutputs = str(len(point_output_times))
    t0 = datetime.datetime.strptime(point_output_times[0],'%Y%m%d %H%M%S')
    t1 = datetime.datetime.strptime(point_output_times[1],'%Y%m%d %H%M%S')
    dt = t1-t0
    output_interval = int(dt.total_seconds())
    if output_interval < 3600:
      output_interval = 3600
    
    # Replace the time information line
    replace_ww3_ounp_inp_line('start date',start.split()[0],start.split()[1],str(output_interval),noutputs)

    for out_type in cfg['out_types']:

      otype = out_types[out_type]
      
      # Replace the output type information lines
      replace_ww3_ounp_inp_line('file prefix' ,otype['prefix'])
      replace_ww3_ounp_inp_line('output type' ,otype['type'])
      replace_ww3_ounp_inp_line('sub-type',otype['subtype'])
  
      # Run the ww3_ounp program
      subprocess.call(['srun','-n','4',pwd+'/ww3_ounp'])
  
  # Move file to data directory
  if not os.path.exists(cfg['data_direc']):
    subprocess.call(['mkdir','-p',cfg['data_direc']])
  subprocess.call('mv *.nc '+cfg['data_direc'],shell=True)
