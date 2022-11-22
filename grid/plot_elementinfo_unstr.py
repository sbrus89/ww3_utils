import numpy as np  
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib.ticker as mticker
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import argparse


def read_gmsh(filename):

    #purpose: this function reads a gmsh file and returns node and element information 
    # additional information about the gmsh file format can be found 
    # in many places including here: http://gmsh.info/dev/doc/texinfo/gmsh.pdf
    # The gmsh format is what the WW3 model uses to define unstructured grids

    #input: 
    # filename - name of gmsh file 

    #output: 
    #xy    -  x/y or lon/lat of nodes (lon points are put in -180 to 180 range)
    #depth - depth value at node points  
    #ect   - element connection table 
    #bnd   - list of boundary nodes 

    #Open gmsh file "filename" 
    f = open(filename,'r')
    lines = f.read().splitlines()

    #skip mesh format lines
    lines.pop(0)
    lines.pop(0)
    lines.pop(0)
    #Skip '$Nodes'
    lines.pop(0)
    #Read number of nodes
    nn = int(lines.pop(0))

    # Initialize node coordinate and depth arrays
    xy = np.zeros((nn,2),dtype=np.double)
    depth = np.zeros((nn,),dtype=np.double)
  
    # Read coordinates and depths
    for i in range(nn):
      line = lines.pop(0).split()
      j = int(line[0])-1
      xy[j,0] = float(line[1]) 
      if xy[j,0] > 180.0:  
         xy[j,0] = xy[j,0] - 360
      xy[j,1] = float(line[2])
      depth[j] = float(line[3])

    #skip '$EndNodes'
    lines.pop(0)
    #Skip '$Elements'
    lines.pop(0)
    #read number of elemends (boundary + elements) 
    ne = int(lines.pop(0))

    #initialize temporary arrays to read in element info
    ecttemp = np.zeros((ne,3),dtype=np.int32)
    bndtemp = np.zeros((ne,1),dtype=np.int32)
    numbnd=0 

    for i in range(ne): 
      #read each element line 
      line = lines.pop(0).split()
      #specify the type of element, boundary (15) or tiangle (2)
      eltype = int(line[1])
      
      if eltype == 15: 
        #store boundary nodes in temp array 
        bndtemp[numbnd] = int(line[5])-1 
        numbnd=numbnd+1
      else:
        #store element connection table in tempoarary array 
        j = int(line[4])-1
        ecttemp[j,0] = int(line[6])-1
        ecttemp[j,1] = int(line[7])-1
        ecttemp[j,2] = int(line[8])-1

    #put element and boundary information in arrays from temp arrays 
    ect = ecttemp[0:j,:]
    bnd = bndtemp[0:numbnd-1] 

    return xy,depth,ect,bnd

def calc_elm_size(xy,ect): 
   
   #purpose: Calculate element size, by calculating the distance between each node on the triangle, then 
   # saving the minimum or maximum distance for each node 

   #input: 
   # xy  -  x/y or lon/lat of nodes 
   # ect - element connection table 

   #output: 
   # distmin(number of nodes)  - the minimum distance between this node and any connected node point
   # distmax(number of nodes)  - the maximum distance between this node and any connected node point


   nn = len(xy)  #number of nodes
   ne = len(ect) #number of elements

   radiusofearth=6378.137 #radius of earth at equator 

   #initialize min and max values
   distmin = np.ones((nn),dtype=np.double)*radiusofearth
   distmax = np.zeros((nn),dtype=np.double)

   for elem in range(ne): 

     #nodes for each element
     i=ect[elem,0]
     j=ect[elem,1] 
     k=ect[elem,2]   
     d2r=np.pi/180 #degrees to radian 
     #calculate the distance between the three points using the haversine formula
     # reference: https://en.wikipedia.org/wiki/Haversine_formula
     distij = 2*radiusofearth*np.arcsin(np.sqrt(np.sin((xy[j,1]-xy[i,1])*d2r/2)**2+np.cos(xy[i,1]*d2r)*np.cos(xy[j,1]*d2r)*np.sin((xy[j,0]-xy[i,0])*d2r/2)**2))
     distjk = 2*radiusofearth*np.arcsin(np.sqrt(np.sin((xy[k,1]-xy[j,1])*d2r/2)**2+np.cos(xy[j,1]*d2r)*np.cos(xy[k,1]*d2r)*np.sin((xy[k,0]-xy[j,0])*d2r/2)**2))
     distki = 2*radiusofearth*np.arcsin(np.sqrt(np.sin((xy[i,1]-xy[k,1])*d2r/2)**2+np.cos(xy[k,1]*d2r)*np.cos(xy[i,1]*d2r)*np.sin((xy[i,0]-xy[k,0])*d2r/2)**2))

     #fill in min/max for each node      
     distmin[i]=np.min([distmin[i],distij,distki]) 
     distmax[i]=np.max([distmax[i],distij,distki])
     distmin[j]=np.min([distmin[j],distij,distjk])
     distmax[j]=np.max([distmax[j],distij,distjk])
     distmin[k]=np.min([distmin[k],distjk,distki])
     distmax[k]=np.max([distmax[k],distjk,distki])

   return distmin,distmax   

def plot_eleminfo(plotdescriptor, xy,ect,distmin,distmax,depth):
 

  print('Grid') 
  print(plotdescriptor) 
  print('Min/max of distmin') 
  print([np.min(distmin),np.max(distmin)])
  print('Min/max of distmax')
  print([np.min(distmax),np.max(distmax)]) 
  print('Min/max of bathy')
  print([np.min(depth),np.max(depth)])

  #Create mask to mask out any element that is going to have ~-180 and 180 nodes 
  # Mask is 0 to include or 1 to exclude 
  # This is a limitation of this plotting routine.  Ideally we'd create new elements to "wrap" values. 

  ne = len(ect) #number of elements
  mask = np.ones((ne,),dtype=np.int32) 
  for elem in range(ne):

     #nodes for each element
     i=ect[elem,0]
     j=ect[elem,1]
     k=ect[elem,2]
     a = [xy[i,0],xy[j,0],xy[k,0]]  
     if len(np.unique(np.sign(a)))==1: 
       mask[elem,] = 0
     else:
       if np.abs(xy[i,0])<10: 
         mask[elem,] = 0
  
  #create one tringulation and use it for multiple plots: 
  triang=tri.Triangulation(xy[:,0],xy[:,1],triangles=ect, mask=mask)
  vpltmin=1
  vpltmax=50

  #make a plot with elements drawn
  fig = plt.figure(figsize=[18.0,9.0])
  ax = fig.add_subplot(1,1,1, projection=ccrs.PlateCarree())
  ax.set_extent([-180,180,-90.0,90.0], crs=ccrs.PlateCarree())
  cf = ax.triplot(triang,'k-', linewidth=0.50)
  ax.add_feature(cfeature.COASTLINE, zorder=101)
  gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=2, color='gray', alpha=0.5, linestyle='--')
  gl.top_labels = False
  gl.right_labels = False
  gl.xlines = False
  gl.ylines = False
  gl.xlocator = mticker.FixedLocator([-180, -120, -60, 0, 60, 120, 180])
  gl.ylocator = mticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])
  gl.xformatter = LongitudeFormatter()
  gl.yformatter = LatitudeFormatter()
  plotname = 'elm_%s.png' % (plotdescriptor)
  plt.savefig(plotname)
  plt.close()

  #plot maximum size of elements 
  fig = plt.figure(figsize=[18.0,9.0])
  ax = fig.add_subplot(1,1,1, projection=ccrs.PlateCarree())
  ax.set_extent([-180,180,-90.0,90.0], crs=ccrs.PlateCarree())
  #var01=ax.tricontourf(triang,distmax,transform=ccrs.PlateCarree())
  var01=ax.tripcolor(triang,distmax,shading='gouraud',vmin=vpltmin, vmax=vpltmax)
  gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=2, color='gray', alpha=0.5, linestyle='--')
  gl.top_labels = False
  gl.right_labels = False
  gl.xlines = False
  gl.ylines = False
  gl.xlocator = mticker.FixedLocator([-180, -120, -60, 0, 60, 120, 180])
  gl.ylocator = mticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])
  gl.xformatter = LongitudeFormatter()
  gl.yformatter = LatitudeFormatter()
  ax.add_feature(cfeature.COASTLINE, zorder=101)
  plt.colorbar(mappable=var01,label="Element size in km")
  plotname = 'maxsize_%s.png' % (plotdescriptor)
  plt.savefig(plotname)
  plt.close()

  #plot minimum size of elements 
  fig = plt.figure(figsize=[18.0,9.0])
  ax = fig.add_subplot(1,1,1, projection=ccrs.PlateCarree())
  ax.set_extent([-180,180,-90.0,90.0], crs=ccrs.PlateCarree())
  var02=ax.tripcolor(triang,distmin,shading='gouraud',vmin=vpltmin, vmax=vpltmax)
  gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=2, color='gray', alpha=0.5, linestyle='--')
  gl.top_labels = False
  gl.right_labels = False
  gl.xlines = False
  gl.ylines = False
  gl.xlocator = mticker.FixedLocator([-180, -120, -60, 0, 60, 120, 180])
  gl.ylocator = mticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])
  gl.xformatter = LongitudeFormatter()
  gl.yformatter = LatitudeFormatter()
  ax.add_feature(cfeature.COASTLINE, zorder=101)
  plt.colorbar(mappable=var02, label="Element size in km")
  plotname = 'minsize_%s.png' % (plotdescriptor)
  plt.savefig(plotname)
  plt.close()


  #plot bathymetry 
  fig = plt.figure(figsize=[18.0,9.0])
  ax = fig.add_subplot(1,1,1, projection=ccrs.PlateCarree())
  ax.set_extent([-180,180,-90.0,90.0], crs=ccrs.PlateCarree())
  var03 = ax.tripcolor(triang,depth) #shading='gouraud',cmap='ocean')
  ax.add_feature(cfeature.COASTLINE, zorder=101)
  gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=2, color='gray', alpha=0.5, linestyle='--')
  gl.top_labels = False
  gl.right_labels = False
  gl.xlines = False
  gl.ylines = False
  gl.xlocator = mticker.FixedLocator([-180, -120, -60, 0, 60, 120, 180])
  gl.ylocator = mticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])
  gl.xformatter = LongitudeFormatter()
  gl.yformatter = LatitudeFormatter()
  plt.colorbar(mappable=var03, label="Bathymetry in m")
  plotname = 'bathy_%s.png' % (plotdescriptor)
  plt.savefig(plotname)
  plt.close()





if __name__ == "__main__":

  ap = argparse.ArgumentParser()
  ap.add_argument('-g', '--gmesh', help="path to gmesh input file", required=True)
  ap.add_argument('-d', '--descriptor', help="descriptor to use for naming plots when saving figures", required=True)
  MyArgs = ap.parse_args()

  xy,depth,ect,bnd = read_gmsh(MyArgs.gmesh)
  distmin,distmax = calc_elm_size(xy,ect)
  plot_eleminfo(MyArgs.descriptor, xy,ect,distmin,distmax,depth)

