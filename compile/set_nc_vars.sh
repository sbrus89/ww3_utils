#!/bin/bash 

export WWATCH3_NETCDF="NC4"
HOST=`hostname -s`
if [[ $HOST = gr* ]] ; then
  export NETCDF_CONFIG="/usr/projects/climate/SHARED_CLIMATE/software/grizzly/netcdf/4.4.1/gcc-5.3.0/bin/nc-config" 
elif [[ $HOST = wf* ]] ; then
  export NETCDF_CONFIG="/usr/projects/climate/SHARED_CLIMATE/software/wolf/netcdf/4.4.0/gcc-4.8.2/bin/nc-config"
fi

echo $WWATCH3_NETCDF
echo $NETCDF_CONFIG


