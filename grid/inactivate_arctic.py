import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import os
import glob

cutoff = 82.0

pwd = os.getcwd()

# Read in mask file
mask_file = glob.glob(pwd+'/*.mask')[0]
mask = np.loadtxt(mask_file)

# Create lon, lat arrays
nx,ny = np.shape(mask)
lon = np.linspace(0.0,360.0,ny)
lat = np.linspace(-90.0,90.0,nx)

# Find indicies of the arctic region
idy, = np.where(lat > cutoff)
idx = np.arange(ny)
index = np.ix_(idy,idx)

# Set points as excluded
mask[index] = 3

np.savetxt(mask_file,mask,fmt='%2.1i',delimiter=' ',newline=' \n')
