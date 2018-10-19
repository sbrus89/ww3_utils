import os
import sys
sys.dont_write_bytecode = True
import yaml
import pprint

def write_ww3_prnc_inp():
  
  pwd = os.getcwd()
  
  f = open(pwd+'/ww3_prnc.config')
  cfg = yaml.load(f)
  pprint.pprint(cfg)

  f = open(pwd+'/ww3_prnc.inp','w')
  
  f.write('$ -------------------------------------------------------------------- $\n')
  f.write('$ WAVEWATCH III Field preprocessor input file                          $\n')
  f.write('$ -------------------------------------------------------------------- $\n')
  f.write('$\n')
  f.write('   '+'  '.join([cfg["field_type"],cfg["format_type"],cfg["time_flag"],cfg["header_flag"]])+'\n')
  f.write('   '+'  '.join([cfg["lon_dimen"],cfg["lat_dimen"]])+'\n')
  f.write('   '+'  '.join([cfg["u_var"],cfg["v_var"]])+'\n')
  f.write('   '+cfg["filename"]+'\n')
  f.write('$ -------------------------------------------------------------------- $\n')
  f.write('$ End of input file                                                    $\n')
  f.write('$ -------------------------------------------------------------------- $\n')
  
if __name__ == '__main__':
  write_ww3_prnc_inp()
