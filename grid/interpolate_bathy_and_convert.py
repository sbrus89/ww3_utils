import fort14
import interpolate_bathymetry
import jigsawpy
import numpy as np
import yaml
import os

########################################################################
########################################################################

def interpolate_bathy(input_file,output_file):
  mesh = fort14.fort14(input_file,verbose=True)
  mesh.read_grid()
  mesh.depth = interpolate_bathymetry.interpolate_global(mesh.xy[:,0], mesh.xy[:,1], filename='etopo1.nc', ele_var='z')
  mesh.write_grid(output_file)

########################################################################
########################################################################

def write_to_vtk(input_file,output_file):
  mesh = fort14.fort14(input_file,verbose=True)
  mesh.read_grid()
  
  #for i in range(mesh.nn):
  #  if mesh.depth[i] < 10.0:
  #    mesh.depth[i] = 10.0
  #  if not np.isfinite(mesh.depth[i]):
  #    print(i,mesh.xy[i,0],mesh.xy[i,1],mesh.depth[i])
  #    mesh.depth[i] = 10.0
  
  msh = jigsawpy.jigsaw_msh_t()
  
  mesh.xy = mesh.xy*np.pi/180.0
  xy = []
  for i in range(mesh.nn):
    xy.append(((mesh.xy[i,0],mesh.xy[i,1]),0))
  msh.vert2 = np.array(xy,dtype=jigsawpy.jigsaw_msh_t.VERT2_t)
  
  ect = []
  mesh.ect = mesh.ect-1
  for i in range(mesh.ne):
    ect.append(((mesh.ect[i,0],mesh.ect[i,1],mesh.ect[i,2]),0))
  msh.tria3 = np.array(ect,dtype=jigsawpy.jigsaw_msh_t.TRIA3_t)
  
  msh.meshID = 'elilipsoid-mesh'
  
  msh.radii = np.full(3, +6371.,
      dtype=jigsawpy.jigsaw_msh_t.REALS_t)
  
  msh.vert3 = np.zeros(msh.vert2.size,
      dtype=jigsawpy.jigsaw_msh_t.VERT3_t)
  
  msh.vert3["coord"] = jigsawpy.S2toR3(
      msh.radii, msh.vert2["coord"])
  
  msh.value = np.array(mesh.depth,dtype=jigsawpy.jigsaw_msh_t.REALS_t)
  print(np.min(mesh.depth))
  print(np.max(mesh.depth))
  print(np.count_nonzero(np.isnan(mesh.depth)))
  
  msh.vert2 = None
  
  jigsawpy.savevtk(output_file, msh)
  
  mesh.xy = mesh.xy*180.0/np.pi
  f = open(output_file,'a')
  #f.write('POINT_DATA '+str(mesh.nn)+'\n')
  f.write('SCALARS lon float 1\n')
  f.write('LOOKUP_TABLE default\n')
  for i in range(mesh.nn):
    f.write(str(mesh.xy[i,0])+'\n')
  #f.write('POINT_DATA '+str(mesh.nn)+'\n')
  f.write('SCALARS lat float 1\n')
  f.write('LOOKUP_TABLE default\n')
  for i in range(mesh.nn):
    f.write(str(mesh.xy[i,1])+'\n')
  f.close()
 
########################################################################
########################################################################

def write_to_gmsh(input_file,output_file): 
  mesh = fort14.fort14(input_file,verbose=True)
  mesh.read_grid()
  for i in range(mesh.nn):
    if mesh.depth[i] < 2.5:
      mesh.depth[i] = 2.5
  mesh.write_gmsh(output_file)
  #x = mesh.xy[:,0]
  #x[x<0] = x[x<0] + 360.0
  #print(mesh.xy)
  #mesh.write_gmsh(output_file+'0-360')

########################################################################
########################################################################
 
if __name__ == '__main__':

  pwd = os.getcwd()

  f = open(pwd+'interpolate_bathy_and_convert.config')
  cfg = yaml.load(f)
  pprint.pprint(cdf)

  interpolate_bathy(cfg['fort14_in'],cfg['fort14_out']) 
  write_to_vtk(cfg['fort14_out'],cfg['vtk_out'])
  write_to_gmsh(cfg['fort14_out'],cfg['gmsh_out'])

