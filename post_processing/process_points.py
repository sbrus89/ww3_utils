import subprocess
import os
import glob

pwd = os.getcwd()

# Set directories
run_direc = '/lustre/scratch4/turquoise/sbrus/WW3_testing/glo_15m/' 
output_direc = pwd+'/model_output/' 

# Output type options
out_types = [{'type':'2','prefix':'ww3.','subtype':'2'},  # mean wave parameters 
             {'type':'2','prefix':'cfsr.','subtype':'1'}] # depth, current, wind

###################################################################################################
###################################################################################################

def replace_ww3_ounp_inp_line(nline,opt1,opt2=None):

  # Find the requested line (nline) in the ww3_ounp input file
  founp = open(pwd+'/ww3_ounp.inp','r')
  lines = founp.read().splitlines()
  line_count = 0
  for n,line in enumerate(lines):
    if line.strip()[0] != '$':
      line_count = line_count + 1
    if line_count == nline:
      break
  
  # Replace the line with the new information (opt1 and opt2) 
  line_info = lines[n].split()
  line_info[0] = opt1
  if opt2:
    line_info[1] = opt2
  lines[n] = '   '+'  '.join(line_info)
  founp.close()
  
  # Re-write the ww3_ounp input file
  founp = open(pwd+'/ww3_ounp.inp','w')
  founp.write('\n'.join(lines))
  founp.close()

###################################################################################################
###################################################################################################

if __name__ == '__main__':

  # Check if the ww3_ounp input file exists
  if not os.path.isfile(pwd+'/ww3_ounp.inp'):
    print 'ww3_ounp.inp not found'
    raise SystemExit(0)
  
  # Link the mod_def.ww3 file to the current directory
  subprocess.call(['ln','-sf',run_direc+'mod_def.ww3',pwd])
  
  # Loop over all out_pnt.ww3.YYYYMMDD_HHMMSS-YYMMDD_HHMMSS files
  pnt_files = sorted(glob.glob(output_direc+'out_pnt.ww3.*'))
  for f in pnt_files:
  
    # Link the out_pnt.ww3 file to the current directory
    subprocess.call(['ln','-sf',f,pwd+'/out_pnt.ww3'])
  
    # Find the start and end dates from the filename
    date_range = f.split(".")[-1]  
    start_date_time = date_range.split('-')[0]
    start_date = start_date_time.split('_')[0]
    start_time = start_date_time.split('_')[1]
    print start_date,start_time
  
    # Replace the time information line
    replace_ww3_ounp_inp_line(1,start_date,start_time)

    for otype in out_types:
      
      # Replace the output type information lines
      replace_ww3_ounp_inp_line(3 ,otype['prefix'])
      replace_ww3_ounp_inp_line(7 ,otype['type'])
      replace_ww3_ounp_inp_line(10,otype['subtype'])
  
      # Run the ww3_ounp program
      subprocess.call(['srun','-n','4',pwd+'/ww3_ounp'])
  
