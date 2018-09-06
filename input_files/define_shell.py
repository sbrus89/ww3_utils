# Input fields to be used
water_lev_flag = 'F'      # These variables are 
current_flag   = 'F'      # assumed to not 
wind_flag      = 'T'      # be homogenous
ice_flag       = 'F'

# Time frame of calculation in yyyymmdd hhmmss format
start_time     = '20050601 000000'
end_time       = '20050630 230000'

# Output server mode
#   0 - No data server process, direct access output from each process 
#   1 - No data server process, output performed by process that performs computations
#   2 - Last process is reserved for all output, does no computing
#   3 - Multiple dedicated output precesses
iostyp         = '1'

#--------------------
# Define output data
#--------------------
# Start and end times in yyyymmdd hhmmss format
# an output interval of 0 turns off that type of output

# Field of mean wave parameters
field_output_start = start_time
field_output_intvl = str(3*3600)
field_output_end   = end_time    
fields = ['HS','TP','DP']        # Assumes namelist selection of fields

# Point output
point_output_start = start_time
point_output_intvl = '3600'      
point_output_end   = end_time
station_file = './stations.txt'

track_output_start = start_time
track_output_intvl = '0'
track_output_end   = end_time

resrt_output_start = start_time
resrt_output_intvl = '0'
resrt_output_end   = end_time

bound_output_start = start_time
bound_output_intvl = '0'
bound_output_end   = end_time

sepfd_output_start = start_time
sepfd_output_intvl = '0'
sepfd_output_end   = end_time

