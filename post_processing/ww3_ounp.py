import os
import yaml
import pprint


def write_ww3_ounp_inp():
  
  pwd = os.getcwd()

  f = open(pwd+'/ww3_ounp.config')
  cfg = yaml.load(f)
  pprint.pprint(cfg)
  
  f = open(pwd+'/ww3_ounp.inp','w')
  
  f.write('$ -------------------------------------------------------------------- $\n')
  f.write('$ WAVEWATCH III NETCDF Point output post-processing                    $\n')
  f.write('$ -------------------------------------------------------------------- $\n')
  f.write('$\n')
  f.write('$ start date, increment, number of outputs\n')
  f.write('   '+'  '.join([cfg['first_output_time'],cfg['increment'],cfg['number_of_outputs']])+'\n')
  f.write('$\n')
  f.write('$ Points requested --------------------------------------------------- $\n')
  for point in cfg['points_requested']:
    f.write('   '+str(point)+'\n')
  f.write('$ mandatory end of list\n')
  f.write('-1\n')
  f.write('$\n')
  f.write('$ -------------------------------------------------------------------- $\n')
  f.write('$ file prefix\n')
  f.write('   '+cfg['file_prefix']+'\n')
  f.write('$ number of characters in date\n')
  f.write('   '+cfg['date_length']+'\n')
  f.write('$ NetCDF version\n')
  f.write('   '+cfg['nc_version']+'\n')
  f.write('$ points in same file, max number of points processed in one pass\n')
  f.write('   '+'  '.join([cfg['points_in_same_file'],cfg['max_points_processed']])+'\n')
  f.write('$ output type\n')
  f.write('   '+cfg['output_type']+'\n')
  f.write('$ flag for global attributes\n')
  f.write('   '+cfg['global_attributes']+'\n')
  f.write('$ flag for dimensions order\n')
  f.write('   '+cfg['dimensions_order']+'\n')
  f.write('$ sub-type\n')
  if cfg['output_type'] == '2':
    f.write('   '+cfg['sub_type']+'\n')
  elif cfg['output_type'] == '1':
    f.write('   '+cfg['sub_type']+' -1 1\n') # Turn off 1-D spectrum
  f.write('$ -------------------------------------------------------------------- $\n')
  f.write('$ End of input file                                                    $\n')
  f.write('$ -------------------------------------------------------------------- $\n')
  
  f.close()
  
if __name__ == '__main__':
  write_ww3_ounp_inp()


