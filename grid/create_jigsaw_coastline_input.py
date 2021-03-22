import os
import numpy as np
import shapefile
import matplotlib.pyplot as plt
plt.switch_backend('agg')



def create_coastline_geometry(shpfiles,outfile,plot_boundaries=False):

  # Initialize 
  pt_list = []
  pt_connect = []
  bou_id = 0
  npt = 0
  
  for shpfile in shpfiles:

  ############################################################
  # Get coastline coordinates and connectivity from shapefile
  ############################################################
  
      sf = shapefile.Reader(shpfile)
      shapes = sf.shapes()
      records = sf.records()
      n = len(sf.shapes())
      for i in range(n):
      
          points = shapes[i].points
          lenpts = len(points)
          if lenpts == 0:
            continue
          #print(records[i])
          #print(lenpts)
      
          #if records[i][0] != '81' and records[i][0] != '214':
          #  continue
      
          idx = slice(0,lenpts-1)
          npt_start = npt
          for pt in points[idx]:
            pt_list.append([np.radians(pt[0]),np.radians(pt[1]),bou_id])
            pt_connect.append([npt, npt+1,bou_id])
            npt = npt+1
          del pt_connect[-1]
          pt_connect.append([npt-1,npt_start,bou_id])
  
  
  ###############################################
  # Plot boundaries
  ###############################################
  
  if plot_boundaries:
    pt_array = np.asarray(pt_list)
    plt.figure()
    for ed in pt_connect:
      x = np.degrees(np.array([pt_array[ed[0],0],pt_array[ed[1],0]]))
      y = np.degrees(np.array([pt_array[ed[0],1],pt_array[ed[1],1]]))
      plt.plot(x,y)
    plt.savefig('outer.png',dpi=500)
    plt.close()
  
  ###############################################
  # Write JIGSAW input 
  ###############################################
  
  f = open(outfile,'w')
  f.write('# Coastline geometry\n')
  f.write('MSHID=3;ELLIPSOID-MESH\n')
  f.write('RADII=6371;6371;6371\n')
  f.write('NDIMS=2\n')
  npt = len(pt_list)
  f.write('POINT='+str(npt)+'\n')
  for pt in pt_list:
    f.write(str(pt[0])+';'+str(pt[1])+';'+str(pt[2])+'\n')
  ned = len(pt_connect)
  f.write('EDGE2='+str(ned)+'\n')
  for ed in pt_connect:
    f.write(str(ed[0])+';'+str(ed[1])+';'+str(ed[2])+'\n')
  f.close()
 

if __name__ == '__main__':

  shpfiles = [
                 "/home/sbrus/run/WW3_unstructured/OceanMesh2D/utilities/GSHHS/l/GSHHS_l_L1.shp",
                 "/home/sbrus/run/WW3_unstructured/OceanMesh2D/utilities/GSHHS/l/GSHHS_l_L6.shp"
             ]

  pwd = os.getcwd()
  create_coastline_geometry(shpfiles,pwd+'/coastlines.msh',plot_boundaries=True)
