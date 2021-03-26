from pathlib import Path
import os
import subprocess
import numpy as np
import xarray
from scipy import interpolate
import jigsawpy
import os
import segment
from create_jigsaw_coastline_input import create_coastline_geometry
from create_jigsaw_initial_input import create_initial_points
import define_hfunction
import calc_distance
from mpas_tools.mesh.creation.jigsaw_to_netcdf import jigsaw_to_netcdf
from mpas_tools.io import write_netcdf
from mpas_tools.mesh.conversion import convert


def ex_8():

    pwd = os.getcwd()
    src_path = pwd 
    dst_path = pwd

    opts = jigsawpy.jigsaw_jig_t()
    
    hfun = jigsawpy.jigsaw_msh_t()
    mesh = jigsawpy.jigsaw_msh_t()
    init = jigsawpy.jigsaw_msh_t()

#------------------------------------ setup files for JIGSAW

    opts.geom_file = \
        str(Path(dst_path)/"geom.msh") # GEOM file
        
    opts.jcfg_file = \
        str(Path(dst_path)/"jcfg.jig") # JCFG file
    
    opts.hfun_file = \
        str(Path(dst_path)/"hfun.msh") # HFUN file

    opts.mesh_file = \
        str(Path(dst_path)/"mesh.msh") # MESH file

    opts.init_file = \
        str(Path(dst_path)/"init.msh") # INIT file
 


#------------------------------------ define JIGSAW geometry
 
    coastlines = False
    if coastlines:

      print('Defining coastline geometry')
      shpfiles = [ 
                     "/home/sbrus/run/WW3_unstructured/OceanMesh2D/utilities/GSHHS/l/GSHHS_l_L1.shp",
                     "/home/sbrus/run/WW3_unstructured/OceanMesh2D/utilities/GSHHS/l/GSHHS_l_L5.shp"
                 ]  

      create_coastline_geometry(shpfiles,opts.geom_file)
    else:
      geom = jigsawpy.jigsaw_msh_t()
      geom.mshID = 'ELLIPSOID-MESH'
      geom.radii = 6371.0*np.ones(3, float)
      jigsawpy.savemsh(opts.geom_file, geom)      
      
    
#------------------------------------ define mesh size function
    
    print('Defining mesh size function')
    ocean_mesh = 'ocean.WC14to60E2r3.200714_scaled.nc'
    km = 1000.0

    lon_min = -180.0
    lon_max = 180.0
    dlon = 0.5
    nlon = int((lon_max-lon_min)/dlon)+1
    lat_min = -90.0
    lat_max = 90.0
    nlat = int((lat_max-lat_min)/dlon)+1
    dlat = 0.5

    
   
    xlon = np.linspace(lon_min,lon_max,nlon)
    ylat = np.linspace(lat_min,lat_max,nlat)
    print(xlon.size,ylat.size)

    hfunction = define_hfunction.cell_widthVsLatLon(xlon,ylat,ocean_mesh)

    hfun.mshID = "ellipsoid-grid"
    hfun.radii = np.full(3, 6371.0, 
        dtype=jigsawpy.jigsaw_msh_t.REALS_t)
    hfun.xgrid = np.radians(xlon)
    hfun.ygrid = np.radians(ylat)
    
    hfun.value = hfunction.astype(jigsawpy.jigsaw_msh_t.REALS_t)

#------------------------------------ specify JIGSAW initial conditions

    print('Specifying initial fixed coordinate locations')
    create_initial_points(ocean_mesh,xlon,ylat,hfunction,opts.init_file)
    jigsawpy.loadmsh(opts.init_file,init)
    jigsawpy.savevtk(str(Path(dst_path)/"init.vtk"), init)

#------------------------------------ set HFUN grad.-limiter

    print('Limiting mesh size function gradients')
    #slope = 0.05
    slope = 0.1
    hfun.slope = np.full(               # |dH/dx| limits
        hfun.value.shape,
        slope, dtype=hfun.REALS_t)

    jigsawpy.savemsh(opts.hfun_file, hfun)

    jigsawpy.cmd.marche(opts, hfun)

#------------------------------------ make mesh using JIGSAW 
    
    opts.hfun_scal = "absolute"
    opts.hfun_hmax = float("inf")   # null HFUN limits
    opts.hfun_hmin = float(+0.00)
    
    opts.mesh_dims = +2             # 2-dim. simplexes

    opts.mesh_eps1 = +.67           # relax edge error
    opts.mesh_rad2 = +1.2
    
    opts.optm_qlim = +9.5E-01       # tighter opt. tol
    opts.optm_iter = +32
    opts.optm_qtol = +1.0E-05
    opts.verbosity = 2

    #jigsawpy.cmd.tetris(opts, 3, mesh)
    jigsawpy.cmd.jigsaw(opts, mesh)
    if coastlines:
      segment.segment(mesh)

#------------------------------------ save mesh for Paraview
    
    print("Writing ex_h.vtk file.")

    jigsawpy.savevtk(
        str(Path(dst_path)/"_hfun.vtk"), hfun)

    jigsawpy.savevtk(
        str(Path(dst_path)/"waves_mesh.vtk"), mesh)

    jigsawpy.savemsh(
        str(Path(dst_path)/"waves_mesh.msh"), mesh)

    jigsaw_to_netcdf(msh_filename='waves_mesh.msh',
                     output_name='waves_mesh_triangles.nc', on_sphere=True,
                     sphere_radius=6371.0)
    write_netcdf(convert(xarray.open_dataset('waves_mesh_triangles.nc'), dir='./',
                         logger=None),
                         'waves_mesh.nc')

    #subprocess.call('./cull_mesh.py --with_critical_passages',shell=True)

    return


if (__name__ == "__main__"): ex_8()



