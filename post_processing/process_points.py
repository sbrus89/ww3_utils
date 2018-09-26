import subprocess
import os
import glob

# Set directories
output_direc = '/lustre/scratch1/turquoise/sbrus/WW3_testing/glo_15m/' 
pwd = os.getcwd()

# Link the mod_def.ww3 file to the current directory
subprocess.call(['ln','-sf',output_direc+'mod_def.ww3',pwd])

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

  # Check if the ww3_ounp input file exists
  if not os.path.isfile(pwd+'/ww3_ounp.inp'):
    print 'ww3_ounp.inp not found'
    raise SystemExit(0)

  # Find the time information line in the ww3_ounp input file
  founp = open(pwd+'/ww3_ounp.inp','r')
  lines = founp.read().splitlines()
  for n,line in enumerate(lines):
    if line.strip()[0] != '$':
      break
  
  # Replace the start date and time to match the out_pnt.ww3 file
  time_info = lines[n].split()
  time_info[0] = start_date
  time_info[1] = start_time
  lines[n] = '   '+'  '.join(time_info)
  founp.close()
  
  # Re-write the ww3_ounp input file
  founp = open(pwd+'/ww3_ounp.inp','w')
  founp.write('\n'.join(lines))
  founp.close()

  # Run the ww3_ounp program
  subprocess.call([pwd+'/ww3_ounp'])

