cmake_minimum_required(VERSION 3.11.3)
project(gmsh2vtk Fortran)
set(sources gmsh2vtk.F90 globals.F90 grid_file_mod.F90 read_gmsh.F90 write_vtk.F90 )
add_executable(gmsh2vtk ${sources})
