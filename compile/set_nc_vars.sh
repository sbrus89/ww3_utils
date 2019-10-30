#!/bin/bash 

export WWATCH3_NETCDF="NC4"
HOST=`hostname -s`
if [[ $HOST = gr* ]] ; then
  export NETCDF_CONFIG="/usr/projects/climate/SHARED_CLIMATE/software/grizzly/netcdf/4.4.1/gcc-5.3.0/bin/nc-config" 
  #export NETCDF_CONFIG="/usr/projects/climate/SHARED_CLIMATE/software/grizzly/netcdf/4.4.1/intel-17.0.1/bin/nc-config" 
elif [[ $HOST = ba* ]] ; then
  export NETCDF_CONFIG="/usr/projects/hpcsoft/toss3/common/netcdf/4.4.0_gcc-7.4.0/bin/nc-config"
elif [[ $HOST = b512 ]] ; then
  export NETCDF_CONFIG="/soft/spack-0.10.0/opt/spack/linux-centos6-x86_64/intel-17.0.1/netcdf-4.4.1.1-mwaopcxxkftrg4jardmcgtw6e7vmxnk4/bin/nc-config"
fi

echo $WWATCH3_NETCDF
echo $NETCDF_CONFIG


