#!/usr/bin/env python
# Phillip J. Wolfram, 10/28/2015
import netCDF4
import numpy as np

def test_assert(logexpr, msg):
    #assert logexpr, msg
    if not logexpr:
        print(msg)
        return True

def test_failed(var, varname, failed):
    failed = failed or test_assert(np.min(var) == 0, '-- ERROR ' + varname + ' field is not bit-zero for min!')
    failed = failed or test_assert(np.max(var) == 0, '-- ERROR ' + varname + ' field is not bit-zero for max!')
    failed = failed or test_assert(np.mean(var) == 0, '-- ERROR ' + varname + ' field is not bit-zero for mean!')
    return failed


def diff(fname1, fname2):
    print('Determining if fields are the same for %s and %s:'%(fname1, fname2))
    print(' ')
    fin1 = netCDF4.Dataset(fname1, 'r')
    fin2 = netCDF4.Dataset(fname2, 'r')
    combinedvars = list(set(fin1.variables) & set(fin2.variables))
    combinedvars.sort()
    print(combinedvars)
    for varname in combinedvars:
        if varname != 'xtime':
          try:
              print('testing %s'%(varname)),
              failed = False
              if varname == 'cellsOnCell' or varname == 'edgesOnCell' or varname == 'verticesOnCell':
                 # for anum, nEdgesOnCell in enumerate(fin1.variables['nEdgesOnCell']):
                 #     var = fin1.variables[varname][anum,:nEdgesOnCell] - fin2.variables[varname][anum,:nEdgesOnCell]
                 #     failed = test_failed(var, varname, failed)
                 #     if failed:
                 #         break
                 pass
              elif varname == 'edgesOnEdge':
                 # for anum, nEdgesOnEdge in enumerate(fin1.variables['nEdgesOnEdge']):
                 #     var = fin1.variables[varname][anum,:nEdgesOnEdge] - fin2.variables[varname][anum,:nEdgesOnEdge]
                 #     if len(var) > 0:
                 #         failed = test_failed(var, varname, failed)
                 #     if failed:
                 #         break
                 pass
              else:
                  shp1 = fin1.variables[varname][:].shape
                  shp2 = fin2.variables[varname][:].shape
                  if shp1 == shp2:
                    var = fin1.variables[varname][:] - fin2.variables[varname][:]
                    failed = test_failed(var, varname, failed)
                  else:
                    var = fin1.variables[varname][-1,:] - fin2.variables[varname][-1,:]
                    failed = test_failed(var, varname, failed)
              if not failed:
                  print('-- ' + varname + ' is bit identical.')
              else:
                #print var
                
                print('min  ', fin1.variables[varname][:].min(), fin2.variables[varname][:].min())
                print('max  ', fin1.variables[varname][:].max(), fin2.variables[varname][:].max())
                print('mean ', fin1.variables[varname][:].mean(), fin2.variables[varname][:].mean())
                print('min  diff ', np.abs(fin1.variables[varname][:] - fin2.variables[varname][:]).min())
                print('max  diff ', np.abs(fin1.variables[varname][:] - fin2.variables[varname][:]).max())
                print('mean diff ', np.abs(fin1.variables[varname][:] - fin2.variables[varname][:]).mean())
                #print 'value', np.vstack((fin1.variables[varname][:], fin2.variables[varname][:]))
          except:
              print('-- WARNING: cannot test ' + varname + ' with values ' + fin1.variables[varname][:] + ' and ' + fin2.variables[varname][:])

if __name__ == "__main__":
    import sys
diff(sys.argv[1], sys.argv[2])
