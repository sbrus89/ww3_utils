import os
import sys
sys.dont_write_bytecode = True
import yaml
import pprint

def write_ww3_shel_inp():
  
  pwd = os.getcwd()
  
  f = open(pwd+'/ww3_shel.config')
  cfg = yaml.load(f)
  pprint.pprint(cfg)

  outputs = ["field_output","point_output","track_output","resrt_output","bound_output","sepfd_output"]

  # Convert output intervals to seconds
  for output in outputs:
    key = output+"_intvl"
    if output == "resrt_output":
      cfg[key] = str(cfg[key]*24*3600)
    else:
      cfg[key] = str(cfg[key]*3600)

  # Set output start/end times to corespond to model start/end times if requested
  for output in outputs:
    for time in ["start","end"]:
      key = output+"_"+time
      if cfg[key] == time+"_time":
        cfg[key] = cfg[time+"_time"]

  f = open(pwd+'/ww3_shel.inp','w')
  
  f.write('$ -------------------------------------------------------------------- $'+'\n')
  f.write('$ WAVEWATCH III shell input file                                       $'+'\n')
  f.write('$ -------------------------------------------------------------------- $'+'\n')
  f.write('$\n')
  f.write('$ Define input to be used ---------------------------------------------$'+'\n')
  f.write('$\n')
  f.write('   '+cfg["water_lev_flag"]+ ' F'+'\n') # These varaiables are
  f.write('   '+cfg["current_flag"]  + ' F'+'\n') # assumed to not 
  f.write('   '+cfg["wind_flag"]     + ' F'+'\n') # be homogenous 
  f.write('   '+cfg["ice_flag"]      + ' F'+'\n')
  f.write('   '+'F F'                 +'\n') # atm momentum 
  f.write('   '+'F F'                 +'\n') # air density
  f.write('   '+'F'                 +'\n') # Assimilation data
  f.write('   '+'F'                 +'\n') # Assimilation data
  f.write('   '+'F'                 +'\n') # Assimilation data
  f.write('$\n')
  f.write('$ Time frame of calculations ------------------------------------------$'+'\n')
  f.write('$\n')
  f.write('   '+cfg["start_time"]+'\n')
  f.write('   '+cfg["end_time"]  +'\n')
  f.write('$\n')
  f.write('$ Define output data --------------------------------------------------$'+'\n')
  f.write('$\n')
  f.write('   '+cfg["iostyp"]+'\n')
  f.write('$\n')
  f.write('$ Fields of mean wave parameters                                        '+'\n')
  f.write('$\n')
  f.write('   '+'  '.join([cfg["field_output_start"],cfg["field_output_intvl"],cfg["field_output_end"]])+'\n')
  if int(cfg["field_output_intvl"]) > 0:
    f.write('   '+'N'+'\n') # Assumes namelist field selection
    f.write('   ')
    for field in cfg["fields"]:
      f.write(field+' ')
    f.write('\n')
  f.write('$\n')
  f.write('$ Point output                                                          '+'\n')
  f.write('$\n')
  f.write('   '+'  '.join([cfg["point_output_start"],cfg["point_output_intvl"],cfg["point_output_end"]])+'\n')
  if int(cfg["point_output_intvl"]) > 0:
    stations = open(cfg["station_file"],'r').read().splitlines()
    for sta in stations:
      f.write('   '+sta+'\n')
    f.write("   0.0  0.0  'STOPSTRING'"+'\n')
  f.write('$\n')
  f.write('$ Output along track                                                    '+'\n')
  f.write('$\n')
  f.write('   '+'  '.join([cfg["track_output_start"],cfg["track_output_intvl"],cfg["track_output_end"]])+'\n')
  f.write('$\n')
  f.write('$ Restart files'+'\n')
  f.write('$\n')
  f.write('   '+'  '.join([cfg["resrt_output_start"],cfg["resrt_output_intvl"],cfg["resrt_output_end"]])+'\n')
  f.write('$\n')
  f.write('$ Boundary data'+'\n')
  f.write('$\n')
  f.write('   '+'  '.join([cfg["bound_output_start"],cfg["bound_output_intvl"],cfg["bound_output_end"]])+'\n')
  f.write('$\n')
  f.write('$ Separated wave field data'+'\n')
  f.write('$\n')
  f.write('   '+'  '.join([cfg["sepfd_output_start"],cfg["sepfd_output_intvl"],cfg["sepfd_output_end"]])+'\n')
  f.write('$\n')
  f.write('$ Homogenous field data ----------------------------------------------$'+'\n')
  f.write('$\n')
  f.write('   '+"'the_end'"+' 0\n')
  f.write('   '+"'STP'"+'\n')
  f.write('$ -------------------------------------------------------------------- $'+'\n')
  f.write('$ End of input file                                                    $'+'\n')
  f.write('$ -------------------------------------------------------------------- $'+'\n')

  f.close()

if __name__ == '__main__':
  write_ww3_shel_inp()
