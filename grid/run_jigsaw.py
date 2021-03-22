from pathlib import Path
import os
import numpy as np
from scipy import interpolate
import jigsawpy
import os
import segment
from create_jigsaw_coastline_input import create_coastline_geometry
from create_jigsaw_initial_input import create_initial_points
import define_hfunction
import calc_distance

def ex_8():

    pwd = os.getcwd()
    src_path = pwd 
    dst_path = pwd

    opts = jigsawpy.jigsaw_jig_t()
    
    geom = jigsawpy.jigsaw_msh_t()
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
 
    opts.verbosity = 2

#------------------------------------ specify JIGSAW initial conditions
    ocean_mesh = 'ocean.WC14to60E2r3.200714_scaled.nc'

    create_initial_points(ocean_mesh,opts.init_file)
    jigsawpy.loadmsh(opts.init_file,init)
    jigsawpy.savevtk(str(Path(dst_path)/"init.vtk"), init)

#------------------------------------ define JIGSAW geometry

    shpfiles = [ 
                   "/home/sbrus/run/WW3_unstructured/OceanMesh2D/utilities/GSHHS/l/GSHHS_l_L1.shp",
                   "/home/sbrus/run/WW3_unstructured/OceanMesh2D/utilities/GSHHS/l/GSHHS_l_L5.shp"
               ]  

    create_coastline_geometry(shpfiles,opts.geom_file)
    
#------------------------------------ define spacing pattern

    km = 1000.0
    hval = 100.0*km

    xlon = np.arange(-180.0, 180.0, 0.5)
    ylat = np.arange(-90.0, 90.0, 0.5)
    print(xlon.size,ylat.size)

    hfunction = define_hfunction.cell_widthVsLatLon(xlon,ylat,ocean_mesh)

    hfun.mshID = "ellipsoid-grid"
    hfun.radii = np.full(3, 6371., 
        dtype=jigsawpy.jigsaw_msh_t.REALS_t)
    hfun.xgrid = np.radians(xlon)
    hfun.ygrid = np.radians(ylat)
    
    hfun.value = hfunction.astype(jigsawpy.jigsaw_msh_t.REALS_t)

#------------------------------------ set HFUN grad.-limiter

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

    #jigsawpy.cmd.tetris(opts, 3, mesh)
    jigsawpy.cmd.jigsaw(opts, mesh)
    segment.segment(mesh)

#------------------------------------ save mesh for Paraview
    
    print("Writing ex_h.vtk file.")

    jigsawpy.savevtk(
        str(Path(dst_path)/"_hfun.vtk"), hfun)

    jigsawpy.savevtk(
        str(Path(dst_path)/"latlon.vtk"), mesh)

    jigsawpy.savemsh(
        str(Path(dst_path)/"latlon.msh"), mesh)

    return


if (__name__ == "__main__"): ex_8()



