water_lev_flag = 'F'      # These variables are 
current_flag   = 'F'      # assumed to not 
wind_flag      = 'T'      # be homogenous
ice_flag       = 'F'

start_time     = '20050601 000000'
end_time       = '20050630 230000'

iostyp         = '1'

field_output_start = start_time
field_output_intvl = str(3*3600)
field_output_end   = end_time

fields = ['HS','TP','DP']

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

