from define_initial import *
import os

pwd = os.getcwd()

f = open(pwd+'/ww3_strt.inp','w')

f.write('$ -------------------------------------------------------------------- $\n')
f.write('$ WAVEWATCH III Initial conditions input file                          $\n')
f.write('$ -------------------------------------------------------------------- $\n')
f.write('$\n')
f.write('   '+itype + '\n')
f.write('$\n')
if itype == 1:
  # Guassian iin frequency and space, cos type in direction
  f.write('   '+'  '.join([]) + '\n')

elif itype == 2:
  # JONSWAP spectrum with Hasselmann et al. (1980) direct. dirstribution
  f.write('   '+'  '.join([]) + '\n')

elif itype == 3:
  # Fetch-limited JONSWAP
  pass

elif itype == 4:
  # User-defined frequency
  f.write('   '+'  '.join([]) + '\n')

elif itype == 5:
  # Start from calm conditions
  pass

f.write('$ -------------------------------------------------------------------- $\n')
f.write('$ End of input file                                                    $\n')
f.write('$ -------------------------------------------------------------------- $\n')

f.close()
