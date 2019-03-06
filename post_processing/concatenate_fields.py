import glob
import subprocess
import datetime
import pprint

############################################################################
############################################################################

# Start and end date range of file to concatenate
start = '08/23/2005' 
end = '08/31/2005'

# Location of files
output_direc = '/Users/sbrus/Projects/DR/glo_15m/fields/august/model_data/'

############################################################################
############################################################################

# Get sorted list of files
files = sorted(glob.glob(output_direc+'*.nc'))
pprint.pprint(files)

# Get datetime object for date range
start_datetime = datetime.datetime.strptime(start,'%m/%d/%Y')
end_datetime = datetime.datetime.strptime(end,'%m/%d/%Y')

# Build ncrcat command
command_list = ['ncrcat']
for f in files:

  # Get datetime object from filename
  date_str = f.split('.')[1]
  date_str = date_str.split('T')[0]
  file_datetime = datetime.datetime.strptime(date_str,'%Y%m%d')
  
  # Add filename to list if it is between date range
  if file_datetime >= start_datetime and file_datetime <= end_datetime:
    command_list.append(f)

# Finish building command
command_list.append('out.nc')
pprint.pprint(command_list)

# Run ncrcat command
subprocess.call(command_list)

 

