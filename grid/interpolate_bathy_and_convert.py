import fort14
import interpolate_bathymetry
import jigsawpy
import numpy as np

def interpolate_bathy(input_file,output_file):
  mesh = fort14.fort14(input_file,verbose=True)
  mesh.read_grid()
  mesh.depth = interpolate_bathymetry.interpolate_global(mesh.xy[:,0], mesh.xy[:,1], filename='etopo1.nc', ele_var='z')
  mesh.write_grid(output_file)

def project_to_sphere(input_file,output_file):
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
 
if __name__ == '__main__':
  #interpolate_bathy('Global_DE_2km.14','fort.14_DE_2km') 
  #project_to_sphere('fort.14_DE_2km','sphere_DE_2km.vtk')
  #write_to_gmsh('fort.14_DE_2km','mesh_DE_2km.msh')

  #interpolate_bathy('clean.14','fort.14_50km_clean') 
  #project_to_sphere('fort.14_50km_clean','sphere_clean.vtk')
  #write_to_gmsh('fort.14_50km_clean','mesh_50km_clean.msh')

  #interpolate_bathy('Global_shallow_4000.14','fort.14_shallow_4000') 
  #project_to_sphere('fort.14_shallow_4000','sphere_shallow_4000.vtk')
  #write_to_gmsh('fort.14_shallow_4000','mesh_shallow_4000.msh')

  #interpolate_bathy('Global_shallow_4000_edit.14','fort.14_shallow_4000_edit') 
  #project_to_sphere('fort.14_shallow_4000_edit','sphere_shallow_4000_edit.vtk')
  #write_to_gmsh('fort.14_shallow_4000_edit','mesh_shallow_4000_edit.msh')

  #interpolate_bathy('Global_shallow_4000_edit_mv_nd.14','fort.14_shallow_4000_edit_mv_nd') 
  #project_to_sphere('fort.14_shallow_4000_edit_mv_nd','sphere_shallow_4000_edit_mv_nd.vtk')
  #write_to_gmsh('fort.14_shallow_4000_edit_mv_nd','mesh_shallow_4000_edit_mv_nd.msh')

  interpolate_bathy('Global_shallow_4000_edit_mv_nd_add_el.14','fort.14_shallow_4000_edit_mv_nd_add_el') 
  project_to_sphere('fort.14_shallow_4000_edit_mv_nd_add_el','sphere_shallow_4000_edit_mv_nd_add_el.vtk')
  write_to_gmsh('fort.14_shallow_4000_edit_mv_nd_add_el','mesh_shallow_4000_edit_mv_nd_add_el.msh')

  #interpolate_bathy('Global_shallow_4000_20.14','fort.14_shallow_4000_20') 
  #project_to_sphere('fort.14_shallow_4000_20','sphere_shallow_4000_20.vtk')
  #write_to_gmsh('fort.14_shallow_4000_20','mesh_shallow_4000_20.msh')

  #interpolate_bathy('Global_half_degree.14','fort.14_half_degree') 
  #project_to_sphere('fort.14_half_degree','sphere_half_degree.vtk')
  #write_to_gmsh('fort.14_half_degree','mesh_half_degree.msh')
