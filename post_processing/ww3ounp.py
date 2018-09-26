from define_points import *
import os

def write_ww3_ounp_inp():
  
  pwd = os.getcwd()
  
  f = open(pwd+'/ww3_ounp.inp','w')
  
  f.write('$ -------------------------------------------------------------------- $\n')
  f.write('$ WAVEWATCH III NETCDF Point output post-processing                    $\n')
  f.write('$ -------------------------------------------------------------------- $\n')
  f.write('$\n')
  f.write('   '+'  '.join([first_output_time,increment,number_of_outputs])+'\n')
  f.write('$\n')
  f.write('$ Points requested --------------------------------------------------- $\n')
  for point in points_requested:
    f.write('   '+point+'\n')
  f.write('$ mandatory end of list\n')
  f.write('-1\n')
  f.write('$\n')
  f.write('$ -------------------------------------------------------------------- $\n')
  f.write('$ file prefix\n')
  f.write('   '+file_prefix+'\n')
  f.write('$ number of characters in date\n')
  f.write('   '+date_length+'\n')
  f.write('$ NetCDF version\n')
  f.write('   '+nc_version+'\n')
  f.write('$ points in same file, max number of points processed in one pass\n')
  f.write('   '+'  '.join([points_in_same_file,max_points_processed])+'\n')
  f.write('$ output type\n')
  f.write('   '+output_type+'\n')
  f.write('$ flag for global attributes\n')
  f.write('   '+global_attributes+'\n')
  f.write('$ flag for dimensions order\n')
  f.write('   '+dimensions_order+'\n')
  f.write('$ sub-type\n')
  f.write('   '+sub_type+'\n')
  f.write('$ -------------------------------------------------------------------- $\n')
  f.write('$ End of input file                                                    $\n')
  f.write('$ -------------------------------------------------------------------- $\n')
  
  f.close()
  
if __name__ == '__main__':
  write_ww3_ounp_inp()


