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
        try:
            print('testing %s'%(varname)),
            failed = False
            shp1 = fin1.variables[varname][:].shape
            shp2 = fin2.variables[varname][:].shape
            if shp1 == shp2:
              var = fin1.variables[varname][:] - fin2.variables[varname][:]
              failed = test_failed(var, varname, failed)
            else:
              print('-- ' + varname + ' not the same shape')
              var = fin1.variables[varname][-1,:] - fin2.variables[varname][-1,:]
              failed = test_failed(var, varname, failed)

            if not failed:
                print('-- ' + varname + ' is bit identical.')
            else:
              print('min  ', np.min(fin1.variables[varname][:]), np.min(fin2.variables[varname][:]))
              print('max  ', np.max(fin1.variables[varname][:]), np.max(fin2.variables[varname][:]))
              print('mean ', np.mean(fin1.variables[varname][:]), np.mean(fin2.variables[varname][:]))
              print('min  diff ', np.min(np.abs(fin1.variables[varname][:] - fin2.variables[varname][:])))
              print('max  diff ', np.max(np.abs(fin1.variables[varname][:] - fin2.variables[varname][:])))
              print('mean diff ', np.mean(np.abs(fin1.variables[varname][:] - fin2.variables[varname][:])))
              #print 'value', np.vstack((fin1.variables[varname][:], fin2.variables[varname][:]))
        except:
            print('-- WARNING: cannot test ' + varname )

if __name__ == "__main__":
    import sys
diff(sys.argv[1], sys.argv[2])
