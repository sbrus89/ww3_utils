#!/bin/bash

HOST=`hostname -s`

if [[ $HOST = gr* ]] ; then
  module purge; module use /usr/projects/climate/SHARED_CLIMATE/modulefiles/all/; module load python/anaconda-2.7-climate; module load gcc/5.3.0 openmpi/1.10.5 netcdf/4.4.1 parallel-netcdf/1.5.0 pio/1.7.2; echo "Loading WaveWatchIII modules for grizzly"
elif [[ $HOST = ba* ]] ; then
  module purge; module load gcc/7.4.0; module load openmpi/2.1.2; module load hdf5-serial/1.8.16; module load netcdf-serial/4.4.0; echo "Loading WaveWatch III modules for badger"
fi

