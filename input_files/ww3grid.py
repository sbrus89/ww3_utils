from define_grid import *
import os

pwd = os.getcwd()

f = open(pwd+'/ww3_grd.inp','w')

f.write('$ -------------------------------------------------------------------- $\n')
f.write('$ WAVEWATCH III Grid preprocessor input file                           $\n')
f.write('$ -------------------------------------------------------------------- $\n')
f.write('$\n')
f.write('  '+grid_name + '\n')
f.write('$\n')
f.write('   '+'  '.join([freq_inc, first_freq, num_freq, num_dir, rel_dir_off]) + '\n')
f.write('$\n')
f.write('   '+' '.join([fldry, flcx, flcy, flcth, flck, flsou]) + '\n')
f.write('$\n')
f.write('   '+'  '.join([max_gl_dt, max_geo_dt, max_spec_dt, min_src_dt])+ '\n')
f.write('$\n')
for ls in namelist:
  f.write('&'+ls+' ')
  for i,var in enumerate(namelist[ls]):
    f.write('   '+var[0]+' = '+var[1])
    if i == len(namelist[ls])-1:
      f.write(' /\n')
    else:
      f.write(', ')
f.write('END OF NAMELISTS\n')
f.write('$\n')
f.write('   '+' '.join([gstrg, flagll, cstrg]) + '\n')
f.write('   '+'  '.join([nx, ny]) + '\n')
f.write('   '+'  '.join([sx, sy, grd_sca_fac]) + '\n')
f.write('   '+'  '.join([x11, y11, crd_sca_fac]) + '\n')
f.write('$ Bathymetry\n')
f.write('   '+'  '.join([lim_bot_dep, min_dep, dep_unit, dep_sca_fac, dep_idla, dep_idfm, dep_frmt, dep_from, dep_name]) + '\n')
f.write('$ Sub-grid information\n')
f.write('   '+'  '.join([obs_unit, obs_sca_fac, obs_idla, obs_idfm, obs_frmt, obs_from, obs_name]) + '\n')
f.write('$ Mask information\n')
f.write('   '+'  '.join([msk_unit, msk_idla, msk_idfm, msk_frmt, msk_from, msk_name]) + '\n')
f.write('$\n')
f.write('   '+4*'0.   '+'0 \n')
f.write('$ -------------------------------------------------------------------- $\n')
f.write('$ End of input file                                                    $\n')
f.write('$ -------------------------------------------------------------------- $\n')

f.close()


