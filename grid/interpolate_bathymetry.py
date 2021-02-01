import jigsawpy
import numpy as np
from scipy import interpolate
from scipy.spatial import Delaunay 
from scipy.spatial import KDTree 
import netCDF4 as nc4
import timeit
import os
import sys
import matplotlib.pyplot as plt
plt.switch_backend('agg')


def read_msh(msh_file):

    mesh = jigsawpy.jigsaw_msh_t()
    jigsawpy.loadmsh(msh_file,mesh)
    nn = len(mesh.point)
    lon_mesh = np.zeros(nn)
    lat_mesh = np.zeros(nn)
    for i in range(nn):
      lon_mesh[i] = mesh.point[i][0][0]
      lat_mesh[i] = mesh.point[i][0][1]

    ne = len(mesh.tria3)
    ect = np.zeros((ne,3))
    for i in range(ne):
      ect[i,:] = mesh.tria3['index'][i][:] 


    return lon_mesh,lat_mesh

########################################################################
########################################################################

def interpolate_global(lon_pts, lat_pts, filename='earth_relief_15s.nc', lon_var='lon', lat_var='lat', ele_var='z'):


    # Open NetCDF data file and read cooordintes
    nc_data = nc4.Dataset(filename, "r")
    lon_data = nc_data.variables[lon_var][:]
    lat_data = nc_data.variables[lat_var][:]

    if np.max(lon_data) > 180.0:
      idx_lon = np.where(lon_pts < 0.0)
      lon_pts[idx_lon] = lon_pts[idx_lon] +  360.0

    # Avoid error with GEBCO lon bounds
    #lon_data[0] = -180.0
    #lon_data[-1] = 180.0

    # Setup interpolation boxes (for large bathymetry datasets)
    n = 2 
    xbox = np.linspace(np.amin(lon_data), np.amax(lon_data), n)
    ybox = np.linspace(np.amin(lat_data), np.amax(lat_data), n)

    dx = xbox[1] - xbox[0]
    dy = ybox[1] - ybox[0]
    boxes = []
    for i in range(n - 1):
        for j in range(n - 1):
            boxes.append(np.asarray(
                [xbox[i], xbox[i + 1], ybox[j], ybox[j + 1]]))

    # Initialize bathymetry
    bathymetry = np.zeros(np.shape(lon_pts))
    #return bathymetry
    bathymetry.fill(np.nan)

    # Interpolate inside each box
    start = timeit.default_timer()
    for i, box in enumerate(boxes):
        print(i + 1, "/", len(boxes))

        # Get points inside box
        lon_idx, = np.where((lon_pts >= box[0]) & (lon_pts <= box[1]))
        lat_idx, = np.where((lat_pts >= box[2]) & (lat_pts <= box[3]))
        idx = np.intersect1d(lon_idx, lat_idx)
        xpts = lon_pts[idx]
        ypts = lat_pts[idx]
        xy_pts = np.vstack((xpts, ypts)).T
        print('  box = ',box)
        print('  points in box: ',xy_pts.size)
        if xy_pts.size == 0: 
          continue

        # Get data inside box (plus a small overlap region)
        overlap = 0.1
        lon_idx, = np.where(
            (lon_data >= box[0] - overlap * dx) & (lon_data <= box[1] + overlap * dx))
        lat_idx, = np.where(
            (lat_data >= box[2] - overlap * dy) & (lat_data <= box[3] + overlap * dy))
        xdata = lon_data[lon_idx]
        ydata = lat_data[lat_idx]
        zdata = nc_data.variables[ele_var][lat_idx, lon_idx]
        testr,testc = np.where(zdata > -30000.0)
        fill_values = zdata.size -  testr.size 
        print('  fill values in data:',fill_values)

        # Interpolate bathymetry onto points
        if fill_values > 0:
          Xdata,Ydata = np.meshgrid(xdata,ydata)
          xpts = Xdata.ravel()
          ypts = Ydata.ravel()
          zpts = zdata.ravel()
          z_idx, = np.where(zpts > -30000.0)
          xpts = xpts[z_idx]
          ypts = ypts[z_idx]
          zpts = zpts[z_idx]

          bathy = interpolate.LinearNDInterpolator((xpts,ypts),zpts)
        else:
          bathy = interpolate.RegularGridInterpolator(
              (xdata, ydata), zdata.T, bounds_error=True, fill_value=np.nan,method='linear')

        bathy_int = bathy(xy_pts)
        bathymetry[idx] = -1.0*bathy_int

    end = timeit.default_timer()
    print(end - start, " seconds")

    if np.max(lon_data) > 180.0:
      lon_pts[idx_lon] = lon_pts[idx_lon] - 360.0


    return bathymetry

########################################################################
########################################################################

def interpolate_csv(lon_pts,lat_pts,bathymetry,filename):

  data = np.genfromtxt(filename, delimiter=',')

  print(lon_pts.size,lat_pts.size)

  xmax = np.amax(data[:,0])
  xmin = np.amin(data[:,0])
  ymax = np.amax(data[:,1])
  ymin = np.amin(data[:,1])
  print(xmin,xmax,ymin,ymax)


  #lon_idx = np.where((lon_pts > xmin) & (lon_pts < xmax))
  #lat_idx = np.where((lat_pts > ymin) & (lat_pts < ymax))
  #idx = np.intersect1d(lon_idx, lat_idx)

  data_pts = np.vstack((data[:,0],data[:,1])).T
  #data_pts = np.vstack((data[:,1],data[:,0])).T
  print(data_pts)
  print(data_pts.shape)
  hull = Delaunay(data_pts)

  mesh_pts = np.vstack((lon_pts,lat_pts)).T
  print(mesh_pts)
  print(mesh_pts.shape)

  idx = np.where(hull.find_simplex(mesh_pts) >= 0)

  xpts = lon_pts[idx]
  ypts = lat_pts[idx]
  xy_pts = np.vstack((xpts,ypts)).T
  print(xy_pts)
  print(xy_pts.shape)

  #bathy = interpolate.NearestNDInterpolator((data[:,0],data[:,1]),data[:,2])
  bathy = interpolate.LinearNDInterpolator((data[:,0],data[:,1]),data[:,2])
  #bathy = interpolate.LinearNDInterpolator((data[:,1],data[:,0]),data[:,2])
  bathy_int = bathy(xy_pts)
  bathymetry[idx] = bathy_int

  return bathymetry
  
########################################################################
########################################################################

if __name__ == "__main__":

    # Open NetCDF mesh file and read mesh points
    msh_file = sys.argv[1]
    lon_mesh, lat_mesh = read_msh(msh_file)
 
    bathymetry = interpolate_global(lon_mesh, lat_mesh, filename='GEBCO_2019.nc', ele_var='elevation')
    
    print(np.amax(bathymetry))
    print(np.amin(bathymetry))
    f = open('bathy.d','w')
    nn = lon_mesh.size
    f.write(str(nn)+'\n')
    for i in range(nn):
      f.write(str(bathymetry[i])+'\n')
    

