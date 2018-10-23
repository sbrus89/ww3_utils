import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
plt.switch_backend('agg')
import glob

dept_file = glob.glob('*.bot')[0]
mask_file = glob.glob('*.mask')[0]
obst_file = glob.glob('*.obst')[0]

dept = np.loadtxt(dept_file)
mask = np.loadtxt(mask_file)
obst = np.loadtxt(obst_file)

nx,ny = np.shape(dept)
print nx,ny
lon = np.linspace(0.0,360.0,ny)
lat = np.linspace(-90.0,90.0,nx)
m = Basemap(projection='cyl', llcrnrlat=-90.0,  urcrnrlat=90.0,
                              llcrnrlon=0.0,    urcrnrlon=360.0, resolution='l')

print "plotting depth"
plt.figure()
plt.contourf(lon,lat,dept)
m.drawcoastlines()
plt.colorbar()
plt.axis('equal')
plt.savefig('depth.png',bbox_inches='tight')

print "plotting land mask"
plt.figure()
plt.contourf(lon,lat,mask)
m.drawcoastlines()
plt.colorbar()
plt.axis('equal')
plt.savefig('mask.png',bbox_inches='tight')

print "plotting x obstructions"
plt.figure()
plt.contourf(lon,lat,obst[1:nx+1,:])
m.drawcoastlines()
plt.colorbar()
plt.axis('equal')
plt.savefig('obst_x.png',bbox_inches='tight')

print "plotting y obstructions"
plt.figure()
plt.contourf(lon,lat,obst[nx:,:])
m.drawcoastlines()
plt.colorbar()
plt.axis('equal')
plt.savefig('obst_y.png',bbox_inches='tight')
