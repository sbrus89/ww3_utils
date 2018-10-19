import os
import sys
sys.dont_write_bytecode = True
import yaml
import pprint

def write_ww3_strt_inp():
  
  pwd = os.getcwd()

  f = open(pwd+'/ww3_strt.config')
  cfg = yaml.load(f)
  pprint.pprint(cfg)

  if cfg["itype"] in ['1','2','4']:
    print "itype not implemented"
    return
  
  f = open(pwd+'/ww3_strt.inp','w')
  
  f.write('$ -------------------------------------------------------------------- $\n')
  f.write('$ WAVEWATCH III Initial conditions input file                          $\n')
  f.write('$ -------------------------------------------------------------------- $\n')
  f.write('$\n')
  f.write('   '+cfg["itype"] + '\n')
  f.write('$\n')
  if cfg["itype"] == '1':
    # Guassian in frequency and space, cos type in direction
    f.write('   '+'  '.join([]) + '\n')
  
  elif cfg["itype"] == '2':
    # JONSWAP spectrum with Hasselmann et al. (1980) direct. dirstribution
    f.write('   '+'  '.join([]) + '\n')
  
  elif cfg["itype"] == '3':
    # Fetch-limited JONSWAP
    pass
  
  elif cfg["itype"] == '4':
    # User-defined frequency
    f.write('   '+'  '.join([]) + '\n')
  
  elif cfg["itype"] == '5':
    # Start from calm conditions
    pass
  
  f.write('$ -------------------------------------------------------------------- $\n')
  f.write('$ End of input file                                                    $\n')
  f.write('$ -------------------------------------------------------------------- $\n')
  
  f.close()

if __name__ == '__main__':
  write_ww3_strt_inp()
