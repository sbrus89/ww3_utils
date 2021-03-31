import pyflann
from scipy import spatial
import numpy as np
import shapefile
import timeit

km = 1000.0

def distance_to_shapefile_points(shpfiles,lon,lat,sphere_radius,reggrid=False,nn_search='flann'):


    # Get coastline coordinates from shapefile
    pt_list = []
    for shpfile in shpfiles:
  
      sf = shapefile.Reader(shpfile)
      shapes = sf.shapes()
      records = sf.records()
      n = len(sf.shapes())

      for i in range(n):
    
          points = shapes[i].points

          if len(points) < 20:
            continue

          for pt in points:
            pt_list.append([pt[0],pt[1]])

    coast_pts = np.radians(np.array(pt_list))


    # Convert coastline points to x,y,z and create kd-tree
    npts = coast_pts.shape[0]
    coast_pts_xyz = np.zeros((npts,3))
    coast_pts_xyz[:, 0], coast_pts_xyz[:, 1], coast_pts_xyz[:, 2] = lonlat2xyz(coast_pts[:, 0], coast_pts[:, 1], sphere_radius)
    if nn_search == "kdtree":
        tree = spatial.KDTree(coast_pts_xyz)
    elif nn_search == "flann":
        flann = pyflann.FLANN()
        flann.build_index(
            coast_pts_xyz,
            algorithm='kdtree',
            target_precision=1.0,
            random_seed=0)

    # Make sure longitude and latitude are in radians
    if np.amax(lon) < 2.1*np.pi:
      lonr = lon
      latr = lat
    else:
      lonr = np.radians(lon)
      latr = np.radians(lat)


    # Convert  backgound grid coordinates to x,y,z and put in a nx x 3 array for kd-tree query
    if reggrid:
      Lon, Lat = np.meshgrid(lonr, latr)
    else:
      Lon = lonr
      Lat = latr
    X, Y, Z = lonlat2xyz(Lon,Lat,sphere_radius)
    pts = np.vstack([X.ravel(), Y.ravel(), Z.ravel()]).T


    # Find distances of background grid coordinates to the coast
    print("   Finding distance")
    start = timeit.default_timer()
    if nn_search == "kdtree":
        d, idx = tree.query(pts)
    elif nn_search == "flann":
        idx, d = flann.nn_index(pts, checks=2000, random_seed=0)
        d = np.sqrt(d)
    end = timeit.default_timer()
    print("   Done")
    print("   " + str(end - start) + " seconds")

    D = np.reshape(d, Lon.shape)

    return D



def lonlat2xyz(lon,lat,R=6371.0*km):

    x = R*np.multiply(np.cos(lon),np.cos(lat))
    y = R*np.multiply(np.sin(lon),np.cos(lat))
    z = R*np.sin(lat)

    return x,y,z
