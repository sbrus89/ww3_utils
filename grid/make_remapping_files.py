import subprocess
import yaml
import os
import pprint

########################################################################
########################################################################

def make_remapping_files(grid1,grid2,grid1_shortname,grid2_shortname,datestamp,reg1,reg2,map_type,nprocs=1):

 
  if map_type == 'conserve':
    map_abbrev = 'aave'
  elif map_type == 'bilinear':
    map_abbrev = 'blin'
  elif map_type == 'patch':
    map_abbrev = 'patc'
  elif map_type == 'neareststod':
    map_abbrev = 'nstod'
  elif map_type == 'nearsetdtos':
    map_abbrev = 'ndtos'
  else:
    print('map type not recognized')
    raise SystemExit(0)

  flags = ' --ignore_unmapped --ignore_degenerate'
  if reg1:
    flags += ' --src_regional'
  if reg2:
    flags += ' --dst_regional'

  map_name = 'map_'+grid1_shortname+'_TO_'+grid2_shortname+'_'+map_abbrev+'.'+str(datestamp)+'.nc'
  subprocess.call('srun -n '+str(nprocs)+' ESMF_RegridWeightGen --source '+grid1+
                                                              ' --destination '+grid2+
                                                              ' --method '+map_type+ 
                                                              ' --weight '+map_name+
                                                              flags, shell=True)

  flags = ' --ignore_unmapped --ignore_degenerate'
  if reg1:
    flags += ' --dst_regional'
  if reg2:
    flags += ' --src_regional'

  map_name = 'map_'+grid2_shortname+'_TO_'+grid1_shortname+'_'+map_abbrev+'.'+str(datestamp)+'.nc'
  subprocess.call('srun -n '+str(nprocs)+' ESMF_RegridWeightGen --source '+grid2+
                                                              ' --destination '+grid1+
                                                              ' --method '+map_type+ 
                                                              ' --weight '+map_name+
                                                              flags, shell=True)
 
########################################################################
########################################################################

if __name__ == '__main__':

  pwd = os.getcwd()

  f = open(pwd+'/make_remapping_files.config')
  cfg = yaml.load(f,Loader=yaml.Loader)
  pprint.pprint(cfg)

  for map_type in cfg['map_types']:

    make_remapping_files(cfg['grid1'],cfg['grid2'],
                         cfg['grid1_shortname'],cfg['grid2_shortname'],
                         cfg['datestamp'],
                         cfg['reg1'],cfg['reg2'],
                         map_type)
  
