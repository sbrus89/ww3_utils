import os
import sys
sys.dont_write_bytecode = True
import yaml
import pprint

def write_ww3_grid_inp():
  
  pwd = os.getcwd()
 
  filenames = ['ww3_grid.config','ww3_grid.config.struc','ww3_grid.config.unst']
  for fname in filenames:
    if os.path.exists(fname):
      f = open(pwd+'/'+fname) 
      break
  cfg = yaml.load(f)
  pprint.pprint(cfg)
  
  
  f = open(pwd+'/ww3_grid.inp','w')
  
  f.write('$ -------------------------------------------------------------------- $\n')
  f.write('$ WAVEWATCH III Grid preprocessor input file                           $\n')
  f.write('$ -------------------------------------------------------------------- $\n')
  f.write('$\n')

  #-----------
  # Grid name 
  #-----------
  f.write('  '+cfg["grid_name"] + '\n')
  f.write('$\n')
 
  #---------------------------
  # Spectral grid information
  #---------------------------
  f.write('   '+'  '.join([cfg["freq_inc"],
                           cfg["first_freq"], 
                           cfg["num_freq"], 
                           cfg["num_dir"], 
                           cfg["rel_dir_off"]]) + '\n')
  f.write('$\n')

  #-------------
  # Model flags
  #-------------
  f.write('   '+' '.join([cfg["fldry"], 
                          cfg["flcx"], 
                          cfg["flcy"], 
                          cfg["flcth"], 
                          cfg["flck"], 
                          cfg["flsou"]]) + '\n')
  f.write('$\n')
 
  #-----------------------
  # Time step parameters
  #-----------------------
  f.write('   '+'  '.join([cfg["max_gl_dt"], 
                           cfg["max_geo_dt"], 
                           cfg["max_spec_dt"], 
                           cfg["min_src_dt"]])+ '\n')
  f.write('$\n')

  #---------------------
  # Namelist parameters
  #---------------------
  for ls in cfg["namelist"]:
    f.write('&'+ls+'\n')
    for i,var in enumerate(cfg["namelist"][ls]):
      f.write('   '+var[0]+' = '+var[1])
      if i == len(cfg["namelist"][ls])-1:
        f.write(' /\n')
      else:
        f.write(',\n')
  f.write('END OF NAMELISTS\n')
  f.write('$\n')

  #---------------------------
  # Grid definition paramters
  #---------------------------
  f.write('   '+' '.join([cfg["gstrg"], 
                          cfg["flagll"], 
                          cfg["cstrg"]]) + '\n')
  if cfg["gstrg"] != "'UNST'":
    f.write('   '+'  '.join([cfg["nx"], 
                             cfg["ny"]]) + '\n')
    f.write('   '+'  '.join([cfg["sx"], 
                             cfg["sy"], 
                             cfg["grd_sca_fac"]]) + '\n')
    f.write('   '+'  '.join([cfg["x11"], 
                             cfg["y11"], 
                             cfg["crd_sca_fac"]]) + '\n')

  #-------------------
  # Bottom depth file
  #-------------------
  f.write('$ Bathymetry\n')
  f.write('   '+'  '.join([cfg["lim_bot_dep"], 
                           cfg["min_dep"], 
                           cfg["dep_unit"], 
                           cfg["dep_sca_fac"], 
                           cfg["dep_idla"], 
                           cfg["dep_idfm"], 
                           cfg["dep_frmt"], 
                           cfg["dep_from"], 
                           cfg["dep_name"]]) + '\n')

  #------------------
  # Obstruction file
  #------------------
  if cfg["gstrg"] != "'UNST'":
    f.write('$ Sub-grid information\n')
    f.write('   '+'  '.join([cfg["obs_unit"],
                             cfg["obs_sca_fac"], 
                             cfg["obs_idla"], 
                             cfg["obs_idfm"], 
                             cfg["obs_frmt"], 
                             cfg["obs_from"], 
                             cfg["obs_name"]]) + '\n')

  #-----------
  # Mask file
  #-----------
  if cfg["gstrg"] != "'UNST'":
    f.write('$ Mask information\n')
    f.write('   '+'  '.join([cfg["msk_unit"], 
                             cfg["msk_idla"],
                             cfg["msk_idfm"],
                             cfg["msk_frmt"],
                             cfg["msk_from"],
                             cfg["msk_name"]]) + '\n')

  #-----------------------
  # Input boundary points
  #-----------------------
  if cfg["gstrg"] == "'UNST'":
    f.write('$ Input boundary points\n')
    f.write('   '+"10 3 1 '(...)' 'PART' 'mapsta.inp'\n")
    f.write('   '+'0   0   F\n') 
    f.write('   '+'0   0   F\n') 
    f.write('   '+'0   0\n') 


  f.write('$\n')
  f.write('   '+4*'0.   '+'0 \n')
  f.write('$ -------------------------------------------------------------------- $\n')
  f.write('$ End of input file                                                    $\n')
  f.write('$ -------------------------------------------------------------------- $\n')
  
  f.close()
  
if __name__ == '__main__':
  write_ww3_grid_inp()


