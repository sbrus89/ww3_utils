cmake_minimum_required(VERSION 3.11.3)
project(grid_utils Fortran)

set(gmsh2vtk_sources gmsh2vtk.F90 globals.F90 grid_file_mod.F90 read_gmsh.F90 write_vtk.F90 )
add_executable(gmsh2vtk ${gmsh2vtk_sources})

set(clean_mesh_sources globals.F90 grid_file_mod.F90 edge_connectivity_mod.F90 fix_elements.F90 clean_mesh.F90)
add_executable(clean_mesh ${clean_mesh_sources})
