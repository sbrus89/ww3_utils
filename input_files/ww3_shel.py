from ww3_shel_config import *
import os

def write_ww3_shel_inp():
  
  pwd = os.getcwd()
  
  f = open(pwd+'/ww3_shel.inp','w')
  
  f.write('$ -------------------------------------------------------------------- $'+'\n')
  f.write('$ WAVEWATCH III shell input file                                       $'+'\n')
  f.write('$ -------------------------------------------------------------------- $'+'\n')
  f.write('$\n')
  f.write('$ Define input to be used ---------------------------------------------$'+'\n')
  f.write('$\n')
  f.write('   '+water_lev_flag+ ' F'+'\n') # These varaiables are
  f.write('   '+current_flag  + ' F'+'\n') # assumed to not 
  f.write('   '+wind_flag     + ' F'+'\n') # be homogenous 
  f.write('   '+ice_flag            +'\n')
  f.write('   '+'F'                 +'\n') # Assimilation data
  f.write('   '+'F'                 +'\n') # Assimilation data
  f.write('   '+'F'                 +'\n') # Assimilation data
  f.write('$\n')
  f.write('$ Time frame of calculations ------------------------------------------$'+'\n')
  f.write('$\n')
  f.write('   '+start_time+'\n')
  f.write('   '+end_time  +'\n')
  f.write('$\n')
  f.write('$ Define output data --------------------------------------------------$'+'\n')
  f.write('$\n')
  f.write('   '+iostyp+'\n')
  f.write('$\n')
  f.write('$ Fields of mean wave parameters                                        '+'\n')
  f.write('$\n')
  f.write('   '+'  '.join([field_output_start,field_output_intvl,point_output_end])+'\n')
  if int(field_output_intvl) > 0:
    f.write('   '+'N'+'\n') # Assumes namelist field selection
    f.write('   ')
    for field in fields:
      f.write(field+' ')
    f.write('\n')
  f.write('$\n')
  f.write('$ Point output                                                          '+'\n')
  f.write('$\n')
  f.write('   '+'  '.join([point_output_start,point_output_intvl,point_output_end])+'\n')
  if int(point_output_intvl) > 0:
    stations = open(station_file,'r').read().splitlines()
    for sta in stations:
      f.write('   '+sta+'\n')
    f.write("   0.0  0.0  'STOPSTRING'"+'\n')
  f.write('$\n')
  f.write('$ Output along track                                                    '+'\n')
  f.write('$\n')
  f.write('   '+'  '.join([track_output_start,track_output_intvl,track_output_end])+'\n')
  f.write('$\n')
  f.write('$ Restart files'+'\n')
  f.write('$\n')
  f.write('   '+'  '.join([resrt_output_start,resrt_output_intvl,resrt_output_end])+'\n')
  f.write('$\n')
  f.write('$ Boundary data'+'\n')
  f.write('$\n')
  f.write('   '+'  '.join([bound_output_start,bound_output_intvl,bound_output_end])+'\n')
  f.write('$\n')
  f.write('$ Separated wave field data'+'\n')
  f.write('$\n')
  f.write('   '+'  '.join([sepfd_output_start,sepfd_output_intvl,sepfd_output_end])+'\n')
  f.write('$\n')
  f.write('$ Homogenous field data ----------------------------------------------$'+'\n')
  f.write('$\n')
  f.write('   '+"'STP'"+'\n')
  f.write('$ -------------------------------------------------------------------- $'+'\n')
  f.write('$ End of input file                                                    $'+'\n')
  f.write('$ -------------------------------------------------------------------- $'+'\n')

  f.close()

if __name__ == '__main__':
  write_ww3_shel_inp()
