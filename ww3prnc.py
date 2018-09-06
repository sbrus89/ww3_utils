from define_input import *
import os

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


