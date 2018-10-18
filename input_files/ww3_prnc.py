from ww3_prnc_config import *
import os

def write_ww3_prnc_inp():
  
  pwd = os.getcwd()
  
  f = open(pwd+'/ww3_prnc.inp','w')
  
  f.write('$ -------------------------------------------------------------------- $\n')
  f.write('$ WAVEWATCH III Field preprocessor input file                          $\n')
  f.write('$ -------------------------------------------------------------------- $\n')
  f.write('$\n')
  f.write('   '+'  '.join([field_type,format_type,time_flag,header_flag])+'\n')
  f.write('   '+'  '.join([lon_dimen,lat_dimen])+'\n')
  f.write('   '+'  '.join([u_var,v_var])+'\n')
  f.write('   '+filename+'\n')
  f.write('$ -------------------------------------------------------------------- $\n')
  f.write('$ End of input file                                                    $\n')
  f.write('$ -------------------------------------------------------------------- $\n')
  
if __name__ == '__main__':
  write_ww3_prnc_inp()
