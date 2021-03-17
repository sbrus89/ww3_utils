import numpy as np
import shapefile
import matplotlib.pyplot as plt
plt.switch_backend('agg')


shpfiles = [
               "/home/sbrus/run/WW3_unstructured/OceanMesh2D/utilities/GSHHS/l/GSHHS_l_L1.shp",
               "/home/sbrus/run/WW3_unstructured/OceanMesh2D/utilities/GSHHS/l/GSHHS_l_L6.shp"
           ]


# Initialize 
pt_list = []
pt_connect = []
bou_id = 0
npt = 0

for shpfile in shpfiles:

    sf = shapefile.Reader(shpfile)
    shapes = sf.shapes()
    records = sf.records()
    n = len(sf.shapes())
    for i in range(n):
    
        points = shapes[i].points
        lenpts = len(points)
        if lenpts == 0:
          continue
        print(records[i])
        print(lenpts)
    
        #if records[i][0] != '81' and records[i][0] != '214':
        #  continue
    
        idx = slice(0,lenpts-1)
        npt_start = npt
        for pt in points[idx]:
          pt_list.append([pt[0],pt[1],bou_id])
          pt_connect.append([npt, npt+1,bou_id])
          npt = npt+1
        del pt_connect[-1]
        pt_connect.append([npt-1,npt_start,bou_id])


###############################################
# Plot boundaries
###############################################

pt_array = np.asarray(pt_list)
plt.figure()
for ed in pt_connect:
  plt.plot([pt_array[ed[0],0],pt_array[ed[1],0]],[pt_array[ed[0],1],pt_array[ed[1],1]])
plt.savefig('outer.png',dpi=500)
plt.close()

###############################################
# Write JIGSAW input 
###############################################

f = open('out.msh','w')
f.write('# test\n')
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

