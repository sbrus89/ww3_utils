from pyevtk.vtk import VtkFile, VtkPolyData
import os
import glob
import numpy as np
import netCDF4
import yaml
import pprint

vec_vars = ['wnd','uss']

########################################################################
########################################################################

def extract_vtk(data_direc, out_dir, out_prefix, variable_list):


    files = sorted(glob.glob(data_direc+'ww3.*'))

    nc_file = open_netcdf(files[0])
    (vertices, connectivity, offsets, M) = get_mesh_information(nc_file)

    build_field_time_series(files,
                            vertices,
                            connectivity, 
                            offsets,
                            out_dir,
                            out_prefix,
                            variable_list,
                            M)

########################################################################
########################################################################

def build_field_time_series(files, vertices, connectivity, offsets, out_dir, out_prefix, variable_list, M):

    outType = 'float64'

    nPolygons = len(offsets)
    nPoints = len(vertices[0])

    # make output directories
    if not os.path.exists(out_dir):
      os.makedirs(out_dir+'/time_series')

    # start the pvd file
    pvd_file = write_pvd_header(out_dir, out_prefix)
    pvd_file.write('<Collection>\n')

    nTimes = len(files)
    for time_index in range(nTimes):
        print(files[time_index])


        # write the header for the vtp file
        vtp_file_prefix = "time_series/{}.{:d}".format(out_prefix, time_index)
        file_name = '{}/{}.vtp'.format(out_dir, vtp_file_prefix)

        nc_file = open_netcdf(files[time_index])
        time = nc_file.variables['time']
        vtk_file = write_vtp_header(out_dir,
                                    vtp_file_prefix,
                                    vertices,
                                    connectivity,
                                    offsets,
                                    nPoints,
                                    nPolygons,
                                    variable_list,
                                    outType,
                                    time)

        # add fields to vtp file
        for var in variable_list:
          print('    '+var)
          if var in vec_vars:
            u = nc_file.variables['u'+var]
            v = nc_file.variables['v'+var]

            uv = np.vstack((np.zeros(u.shape),u,v)).T

            xyz = np.einsum('ijk,ik->ij',M,uv)
            x = np.squeeze(xyz[:,0])
            y = np.squeeze(xyz[:,1])
            z = np.squeeze(xyz[:,2])

            x[x>1e10] = np.nan
            y[y>1e10] = np.nan
            z[z>1e10] = np.nan

            vtk_file.appendData(x.astype(outType))
            vtk_file.appendData(y.astype(outType))
            vtk_file.appendData(z.astype(outType))
          else:
            field = np.squeeze(nc_file.variables[var])
            field[field > 1e10] = np.nan
            field = field.astype(outType)
            vtk_file.appendData(field)
        vtk_file.save()

        # add time step to pdv file
        pvd_file.write('<DataSet timestep="{:.16f}" group="" '
                       'part="0"\n'.format(time_index))
        pvd_file.write('\tfile="{}.vtp"/>\n'.format(vtp_file_prefix))


    # finish the pdv file
    pvd_file.write('</Collection>\n')
    pvd_file.write('</VTKFile>\n')
    pvd_file.close()  # }}}

########################################################################
########################################################################

def open_netcdf(file_name):
    nc_file = netCDF4.Dataset(file_name, 'r')
    # turn off auto mask (if applicable)
    try:
        nc_file.set_auto_mask(False)
    except AttributeError:
        pass
    return nc_file

########################################################################
########################################################################

def write_pvd_header(path, prefix):  
    pvd_file = open('{}/{}.pvd'.format(path, prefix), 'w')
    pvd_file.write('<?xml version="1.0"?>\n')
    pvd_file.write('<VTKFile type="Collection" version="0.1"\n')
    pvd_file.write('\tbyte_order="LittleEndian"\n')
    pvd_file.write('\tcompressor="vtkZLibDataCompressor">\n')
    return pvd_file  

########################################################################
########################################################################

def write_vtp_header(path, prefix, 
                     vertices, connectivity,
                     offsets, nPoints, nPolygons, variable_list, outType, xtime=None): 
    vtkFile = VtkFile("{}/{}".format(path, prefix), VtkPolyData)

    #if xtime is not None:
    #    vtkFile.openElement(str("metadata"))
    #    vtkFile.openElement(str("xtime"))
    #    vtkFile.xml.addText(str(xtime))
    #    vtkFile.closeElement(str("xtime"))
    #    vtkFile.closeElement(str("metadata"))

    vtkFile.openElement(vtkFile.ftype.name)
    vtkFile.openPiece(npoints=nPoints, npolys=nPolygons)

    vtkFile.openElement(str("Points"))
    vtkFile.addData(str("points"), vertices)
    vtkFile.closeElement(str("Points"))

    vtkFile.openElement(str("Polys"))
    vtkFile.addData(str("connectivity"), connectivity)
    vtkFile.addData(str("offsets"), offsets)
    vtkFile.closeElement(str("Polys"))

    vtkFile.openData(str("Point"),
                     scalars='solution')
    for var in variable_list:
        if var in vec_vars:
          vtkFile.addHeader('x'+var, outType, nPoints, 1)
          vtkFile.addHeader('y'+var, outType, nPoints, 1)
          vtkFile.addHeader('z'+var, outType, nPoints, 1)
        else:
          vtkFile.addHeader(var, outType, nPoints, 1)
    vtkFile.closeData(str("Point"))

    vtkFile.closePiece()
    vtkFile.closeElement(vtkFile.ftype.name)

    vtkFile.appendData(vertices)
    vtkFile.appendData(connectivity)
    vtkFile.appendData(offsets)

    return vtkFile  

########################################################################
########################################################################

def get_mesh_information(nc_file):  

    lon = np.radians(nc_file.variables['longitude'][:])
    lat = np.radians(nc_file.variables['latitude'][:])
    nnode = lon.shape[0]

    R = 6371.0
    X = R*np.cos(lat)*np.cos(lon)
    Y = R*np.cos(lat)*np.sin(lon)
    Z = R*np.sin(lat)
    vertices = (X,Y,Z)

    ect = nc_file.variables['tri'][:, :] - 1

    connectivity = np.array(ect.ravel(), int)
    nelem = ect.shape[0]
    offsets = np.array(3*np.arange(1, nelem+1), int)

    M = Tmat(nnode,lon,lat)

    return vertices, connectivity, offsets, M

########################################################################
########################################################################

def Tmat(nn,lon,lat):

  M = np.zeros((nn,3,3))
  M[:,0,0] = np.cos(lat)*np.cos(lon)
  M[:,1,0] = np.cos(lat)*np.sin(lon)
  M[:,2,0] = np.sin(lat)

  M[:,0,1] = -np.sin(lon)
  M[:,1,1] = np.cos(lon)
  M[:,2,1] = 0.0

  M[:,0,2] = -np.sin(lat)*np.cos(lon)
  M[:,1,2] = -np.sin(lat)*np.sin(lon)
  M[:,2,2] = np.cos(lat)

  return M

########################################################################
########################################################################

if __name__ == '__main__':

    pwd = os.getcwd()
    
    # Read config file  
    f = open(pwd+'/write_vtk.config')
    cfg = yaml.load(f)
    pprint.pprint(cfg)

    extract_vtk(cfg['data_dir'], cfg['out_dir'], cfg['out_prefix'], cfg['variable_list'])

