import subprocess
import os
import glob
import yaml
import pprint
import ww3_log
import datetime

pwd = os.getcwd()

# Output type options
point_out_types = {'wave':{'type':'2','prefix':'ww3.','subtype':'2'},  # mean wave parameters 
                   'met': {'type':'2','prefix':'cfsr.','subtype':'1'}, # depth, current, wind
                   'spec':{'type':'1','prefix':'spec.','subtype':'3'}} 
###################################################################################################
###################################################################################################

def replace_ww3_ounp_inp_line(filename,comment,opt1=None,opt2=None,opt3=None,opt4=None,ls=None):

  # Find the requested line (nline) in the ww3_ounp input file
  founp = open(pwd+'/'+filename,'r')
  lines = founp.read().splitlines()
  for n,line in enumerate(lines):
    if line.find(comment) > 0:
      break
  n = n+1
  
  # Replace the line with the new information (opt1 and opt2) 
  line_info = lines[n].split()
  if ls:
    line_info = ls
  if opt1:
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
  founp = open(pwd+'/'+filename,'w')
  founp.write('\n'.join(lines))
  founp.close()

###################################################################################################
###################################################################################################

def write_ounf_inp():

  contents = ['$ -------------------------------------------------------------------- $',
              '$ WAVEWATCH III NETCDF Grid output post-processing                     $',   
              '$ -------------------------------------------------------------------- $',
              '$                                                                       ',
              '$ start date, increment, number of outputs                              ',
              '   20050601  000000  3600  721                                          ',
              '$                                                                       ',
              '$ Namelist type selection--------------------------------------------- $',
              '    N                                                                   ',
              '$ Fields requested --------------------------------------------------- $',
              '   HS  FP  DP                                                           ',
              '$                                                                       ',
              '$ -------------------------------------------------------------------- $',
              '$ NetCDF version                                                        ',
              '   4 4                                                                  ',
              '$ swell partitions                                                      ',
              '   0  1  2                                                              ',
              '$ variables in same file                                                ',
              '   T                                                                    ',
              '$ file prefix                                                           ',
              '   ww3.                                                                 ',
              '$ number of characters in date                                          ',
              '   10                                                                   ',
              '$ ix and iy ranges                                                      ',
              '   0  1000000  0  1000000                                               ',
              '$ -------------------------------------------------------------------- $',
              '$ End of input file                                                    $',
              '$ -------------------------------------------------------------------- $']
  f = open(pwd+'/ww3_ounf.inp','w')
  f.write('\n'.join(contents))

###################################################################################################
###################################################################################################

def write_ounp_inp():

  contents = ['$ -------------------------------------------------------------------- $',
              '$ WAVEWATCH III NETCDF Point output post-processing                    $',
              '$ -------------------------------------------------------------------- $',
              '$                                                                       ',
              '$ start date, increment, number of outputs                              ',
              '   20050601  000000  3600  721                                          ',
              '$                                                                       ',
              '$ Points requested --------------------------------------------------- $',
              '$ mandatory end of list                                                 ',
              '-1                                                                      ',
              '$                                                                       ',
              '$ -------------------------------------------------------------------- $',
              '$ file prefix                                                           ',
              '   cfsr.                                                                ',
              '$ number of characters in date                                          ',
              '   10                                                                   ',
              '$ NetCDF version                                                        ',
              '   4                                                                    ',
              '$ points in same file, max number of points processed in one pass       ',
              '   T  150                                                               ',
              '$ output type                                                           ',
              '   2                                                                    ',
              '$ flag for global attributes                                            ',
              '   0                                                                    ',
              '$ flag for dimensions order                                             ',
              '   T                                                                    ',
              '$ sub-type                                                              ',
              '   1                                                                    ',
              '$ -------------------------------------------------------------------- $',
              '$ End of input file                                                    $',
              '$ -------------------------------------------------------------------- $']
  f = open(pwd+'/ww3_ounp.inp','w')
  f.write('\n'.join(contents))

###################################################################################################
###################################################################################################

if __name__ == '__main__':

  f = open(pwd+'/process_output.config')
  cfg = yaml.load(f,yaml.Loader)
  pprint.pprint(cfg)

  
  # Link the mod_def.ww3 file to the current directory
  subprocess.call(['ln','-sf',cfg['run_direc']+'mod_def.ww3',pwd])

  for out_type in cfg['out_types']:

    if out_type == 'fields':
      filename = 'out_grd'
      exe = 'ww3_ounf'
    elif out_type == 'points':
      filename = 'out_pnt'
      exe = 'ww3_ounp'
    else:
      print('output type not recongnized')
      raise SystemExit(0)



    # Check if the ww3_ounp input file exists
    if not os.path.isfile(pwd+'/'+exe+'.inp'):
      if exe == 'ww3_ounf':
        write_ounf_inp()
      else:
        write_ounp_inp()
    
    # Loop over all out_pnt.ww3.YYYYMMDD_HHMMSS-YYMMDD_HHMMSS files
    pnt_files = sorted(glob.glob(cfg['output_direc']+filename+'.ww3*'))
    log_files = sorted(glob.glob(cfg['log_direc']+'log.ww3*'))
    for i in range(len(pnt_files)):
    
      f = pnt_files[i]
      # Link the out_pnt.ww3 file to the current directory
      subprocess.call(['ln','-sf',f,pwd+'/'+filename+'.ww3'])
    
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
      replace_ww3_ounp_inp_line(exe+'.inp','start date',opt1=start.split()[0],opt2=start.split()[1],opt3=str(output_interval),opt4=noutputs)

      if out_type == 'points':
  
        for point_out_type in cfg['point_out_types']:
  
          otype = point_out_types[point_out_type]
          
          # Replace the output type information lines
          replace_ww3_ounp_inp_line(exe+'.inp','file prefix' ,opt1=otype['prefix'])
          replace_ww3_ounp_inp_line(exe+'.inp','output type' ,opt1=otype['type'])
          replace_ww3_ounp_inp_line(exe+'.inp','sub-type',opt1=otype['subtype'])
    
          # Run the ww3_ounp program
          subprocess.call(['srun','-n','4',pwd+'/ww3_ounp'])

      elif out_type == 'fields':

        #Run the ww3_ounf program
        replace_ww3_ounp_inp_line(exe+'.inp','Fields requested' ,ls=cfg['fields_requested'])
        subprocess.call([pwd+'/ww3_ounf'])
    
    # Move file to data directory
    if not os.path.exists(cfg['data_direc']+out_type):
      subprocess.call(['mkdir','-p',cfg['data_direc']+out_type])
    subprocess.call('mv *.nc '+cfg['data_direc']+out_type,shell=True)
