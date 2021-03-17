
from pathlib import Path
import os
import numpy as np
from scipy import interpolate
import jigsawpy
import os
import segment

def ex_8():

# DEMO-8: build regional meshes via stereographic projection

    pwd = os.getcwd()
    src_path = pwd 
    dst_path = pwd

    opts = jigsawpy.jigsaw_jig_t()
    
    geom = jigsawpy.jigsaw_msh_t()
    hfun = jigsawpy.jigsaw_msh_t()
    mesh = jigsawpy.jigsaw_msh_t()
    mesh_edit = jigsawpy.jigsaw_msh_t()

#------------------------------------ setup files for JIGSAW

    opts.geom_file = \
        str(Path(dst_path)/"geom.msh") # GEOM file
        
    opts.jcfg_file = \
        str(Path(dst_path)/"jcfg.jig") # JCFG file
    
    opts.hfun_file = \
        str(Path(dst_path)/"hfun.msh") # HFUN file

    opts.mesh_file = \
        str(Path(dst_path)/"mesh.msh") # MESH file

#------------------------------------ define JIGSAW geometry

    jigsawpy.loadmsh(
        str(Path(src_path)/"out.msh"), geom)
    geom.point["coord"][:,:] *= np.pi/180.
    jigsawpy.savemsh(opts.geom_file, geom)

    xmin = -180.0
    xmax = 180.0
    xlon = np.arange(xmin, xmax, 0.5)

    ymin = -90.0
    ymax = 90.0
    ylat = np.arange(ymin, ymax, 0.5)

    print(xlon.size,ylat.size)
    
#------------------------------------ define spacing pattern

    km = 1000.0
    hval = 100.0*km

    hfunction = np.full(
        (ylat.size, xlon.size), hval)

    hfun.mshID = "ellipsoid-grid"
    hfun.radii = np.full(3, 6371., 
        dtype=jigsawpy.jigsaw_msh_t.REALS_t)
    hfun.xgrid = xlon * np.pi / 180.
    hfun.ygrid = ylat * np.pi / 180.
    
    hfun.value = hfunction.astype(jigsawpy.jigsaw_msh_t.REALS_t)/km

#------------------------------------ set HFUN grad.-limiter

    hfun.slope = np.full(               # |dH/dx| limits
        hfun.value.shape,
        +0.050, dtype=hfun.REALS_t)

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



