import os
import subprocess

run_name = '20211221.AWAV_WCYCLSSP585_CMIP6.ne30_ECv3_wQU225EC60to30.anvil'
post_processing_dir = './post_processing'
obs_dir = '/home/sbrus/run/WW3_testing/glo_30m/post_processing/obs_data/'
output_dir = '/lcrc/group/acme/sbrus/scratch/anvil/'+run_name+'/run/'
run_dir = output_dir
log_dir = output_dir
nc_output_dir = '/lcrc/group/acme/sbrus/scratch/anvil/'+run_name+'/post_processing/model_data/'

ww3_utils = './ww3_utils/post_processing/'
ww3_exe_dir = '/gpfs/fs1/home/sbrus/E3SM/components/ww3/src/source/WW3/model/exe/'

########################################################################
########################################################################

def replace_config_file_line(direc,config_file,option,replacement):

  pwd = os.getcwd()
  os.chdir(direc)
  f = open(config_file,'r')
  lines = f.read().splitlines()
  f.close()
  for i,line in enumerate(lines):
    if line.find(option) >= 0:
      lines[i] = option+" : '"+replacement+"'"
  
  f = open(config_file,'w')
  f.write('\n'.join(lines))
  f.close()
  os.chdir(pwd)

########################################################################
########################################################################

commands = [['mkdir','-p',post_processing_dir],
            ['ln', '-srf', ww3_utils+'ww3_ounp.py',            post_processing_dir],
            ['ln', '-srf', ww3_utils+'ww3_ounf.py',            post_processing_dir],
            ['ln', '-srf', ww3_utils+'process_output.py',      post_processing_dir],
            ['ln', '-srf', ww3_utils+'plot_points.py',         post_processing_dir],
            ['ln', '-srf', ww3_utils+'plot_fields.py',         post_processing_dir],
            ['ln', '-srf', ww3_utils+'write_vtk.py',           post_processing_dir],
            ['cp',         ww3_utils+'ww3_ounp.config',        post_processing_dir],
            ['cp',         ww3_utils+'ww3_ounf.config',        post_processing_dir],
            ['cp',         ww3_utils+'process_output.config',  post_processing_dir],
            ['cp',         ww3_utils+'plot_points.config',     post_processing_dir],
            ['cp',         ww3_utils+'plot_fields.config',     post_processing_dir],
            ['cp',         ww3_utils+'write_vtk.config',       post_processing_dir],
            ['cp',         ww3_exe_dir+'ww3_ounp',             post_processing_dir],
            ['cp',         ww3_exe_dir+'ww3_ounf',             post_processing_dir],
            ['ln', '-srf', obs_dir,                            post_processing_dir]]
            

for cmd in commands:
  print(' '.join(cmd))
  subprocess.call(cmd)

replace_config_file_line(post_processing_dir,'process_output.config','output_direc',output_dir)
replace_config_file_line(post_processing_dir,'process_output.config','run_direc',run_dir)
replace_config_file_line(post_processing_dir,'process_output.config','log_direc',log_dir)
replace_config_file_line(post_processing_dir,'process_output.config','data_direc',nc_output_dir+'points/')
os.chdir(post_processing_dir)
subprocess.call(['python', 'ww3_ounf.py'])
subprocess.call(['python', 'ww3_ounp.py'])
