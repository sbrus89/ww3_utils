import os
import yaml
import pprint


def write_ww3_ounf_inp():
  
  pwd = os.getcwd()

  f = open(pwd+'/ww3_ounf.config')
  cfg = yaml.load(f,yaml.Loader)
  pprint.pprint(cfg)
  
  f = open(pwd+'/ww3_ounf.inp','w')
  
  f.write('$ -------------------------------------------------------------------- $\n')
  f.write('$ WAVEWATCH III NETCDF Grid output post-processing                     $\n')
  f.write('$ -------------------------------------------------------------------- $\n')
  f.write('$\n')
  f.write('   '+'  '.join([cfg['first_output_time'],
                           cfg['increment'],
                           cfg['number_of_outputs']])+'\n')
  f.write('$\n')
  f.write('$ Fields requested --------------------------------------------------- $\n')
  f.write('    N\n')                        # Assumes namelist specification
  f.write('   '+'  '.join(cfg['fields_requested'])+'\n')
  f.write('$\n')
  f.write('$ -------------------------------------------------------------------- $\n')
  f.write('$ NetCDF version\n')
  f.write('   '+cfg['nc_version']+' 4\n')   # Assumes REAL variable type
  f.write('$ swell partitions\n')
  f.write('   '+'  '.join(cfg['swell_partitions'])+'\n')
  f.write('$ variables in same file\n')
  f.write('   '+cfg['variables_in_same_file']+'\n')
  f.write('$ file prefix\n')
  f.write('   '+cfg['file_prefix']+'\n')
  f.write('$ number of characters in date\n')
  f.write('   '+cfg['date_length']+'\n')
  f.write('$ ix and iy ranges\n')
  f.write('   '+'  '.join([cfg['ix_range'][0],
                           cfg['ix_range'][1],
                           cfg['iy_range'][0],
                           cfg['iy_range'][1]])+'\n')
  f.write('$ -------------------------------------------------------------------- $\n')
  f.write('$ End of input file                                                    $\n')
  f.write('$ -------------------------------------------------------------------- $\n')
  
  f.close()
  
if __name__ == '__main__':
  write_ww3_ounf_inp()

