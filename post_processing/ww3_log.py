import os

#############################################################################################################

def find_output_times(log_file):

  if not os.path.isfile(log_file):
    print('log.ww3 file does not exist')
    raise SystemExit(0)

  f = open(log_file,'r')
  lines = f.read().splitlines()
 
  output_intervals = False
  start_time = ''
  end_time = ''
  restart_output_times = []
  gridded_output_times = []
  point_output_times = []
  n = 0
  for i,line in enumerate(lines):

    n = n + 1

    # Find lines inside output table
    if line.find('--------|------|---------------------|-------------------|---------------|') >= 0:
      output_intervals = True
      n = 0 
      continue
    if line.find('--------+------+---------------------+-------------------+---------------+') >= 0 and n > 1: # For whatever reason, for E3SM output this line is  
      output_intervals = False                                                                                 # right afer the ----|----|----| line with output 
    #elif n == 1:                                                                                               # following it. This prevents if from being flaged
    #  continue                                                                                                 # as the end of the file    

    # Process lines inside output table
    if output_intervals == True:
      cols = line.split('|')
      print(cols)
      
      # Find date and time
      date_time = cols[2].strip().split()
      if len(date_time) > 1:
        date = date_time[0]
        time = date_time[1]
      else:
        time = date_time[0]
      date_time_formatted = ''.join(date.split('/')) + ' ' + ''.join(time.split(':'))

      # Find start time
      if start_time == '':
        start_time = date_time_formatted

      output = cols[4]

      # Find restart output
      if output[7] == 'X' or output[7] == 'L':
        restart_output_times.append(date_time_formatted)
      # Find grided output
      if output[1] == 'X' or output[1] == 'L':
        gridded_output_times.append(date_time_formatted)
      # Find point output
      if output[3] == 'X' or output[3] == 'L':
        point_output_times.append(date_time_formatted)

  end_time = date_time_formatted

  return restart_output_times,gridded_output_times,point_output_times,start_time,end_time

#############################################################################################################

if __name__ == '__main__':
  pass





