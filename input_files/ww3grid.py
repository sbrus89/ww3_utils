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
f.write('   '+'  '.join([str(freq_inc), str(first_freq), str(num_freq), str(num_dir), str(rel_dir_off)]) + '\n')
f.write('$\n')
f.write('   '+' '.join([fldry, flcx, flcy, flcth, flck, flsou]) + '\n')
f.write('$\n')
f.write('   '+'  '.join([str(max_gl_dt), str(max_geo_dt), str(max_spec_dt), str(min_src_dt)])+ '\n')
f.write('$\n')
for ls in namelist:
  f.write('&'+ls+' ')
  for i,var in enumerate(namelist[ls]):
    f.write('   '+var[0]+' = '+str(var[1]))
    if i == len(namelist[ls])-1:
      f.write(' /\n')
    else:
      f.write(', ')
f.write('END OF NAMELISTS\n')
f.write('$\n')
f.write('   '+' '.join([gstrg, flagll, cstrg]) + '\n')
f.write('   '+'  '.join([str(nx), str(ny)]) + '\n')
f.write('   '+'  '.join([str(sx), str(sy), str(grd_sca_fac)]) + '\n')
f.write('   '+'  '.join([str(x11), str(y11), str(crd_sca_fac)]) + '\n')
f.write('$ Bathymetry\n')
f.write('   '+'  '.join([str(lim_bot_dep), str(min_dep), str(dep_unit), str(dep_sca_fac), str(dep_idla), str(dep_idfm), dep_frmt, dep_from, dep_name]) + '\n')
f.write('$ Sub-grid information\n')
f.write('   '+'  '.join([str(obs_unit), str(obs_sca_fac), str(obs_idla), str(obs_idfm), str(obs_frmt), str(obs_from), str(obs_name)]) + '\n')
f.write('$ Mask information\n')
f.write('   '+'  '.join([str(msk_unit), str(msk_idla), str(msk_idfm), str(msk_frmt), str(msk_from), str(msk_name)]) + '\n')
f.write('$\n')
f.write('   '+4*'0.   '+'0 \n')
f.write('$ -------------------------------------------------------------------- $\n')
f.write('$ End of input file                                                    $\n')
f.write('$ -------------------------------------------------------------------- $\n')

f.close()


