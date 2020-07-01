#!/bin/bash

############################################
# Edit this area
############################################

GRID1=ocean.RRS.15-5km_scrip_150722.nc
GRID2=T62_040121.nc
GRID1_SHORTNAME=oRRS15to5
GRID2_SHORTNAME=T62
DATESTAMP=150722
REG1=0
REG2=0
NPROCS=$1

echo $GRID1

############################################
# Don't edit below here
############################################

MAPTYPES="conserve bilinear patch neareststod nearestdtos"

for TYPE in ${MAPTYPES}
do
	if [ "${TYPE}" == "conserve" ]; then
		MAP_ABBREV="aave"
	elif [ "${TYPE}" == "bilinear" ]; then
		MAP_ABBREV="blin"
	elif [ "${TYPE}" == "patch" ]; then
		MAP_ABBREV="patc"
	elif [ "${TYPE}" == "neareststod" ]; then
		MAP_ABBREV="nstod"
	elif [ "${TYPE}" == "nearestdtos" ]; then
		MAP_ABBREV="ndtos"
	fi

	FLAGS="--ignore_unmapped"

	if [ ${REG1} == 1 ]; then
		FLAGS="${FLAGS} --src_regional"
	fi

	if [ ${REG2} == 1 ]; then
		FLAGS="${FLAGS} --dst_regional"
	fi

	MAP_NAME="map_${GRID1_SHORTNAME}_TO_${GRID2_SHORTNAME}_${MAP_ABBREV}.${DATESTAMP}.nc"
	ESMF_RegridWeightGen --source ${GRID1} --destination ${GRID2} --method ${TYPE} --weight ${MAP_NAME} ${FLAGS} 


	FLAGS="--ignore_unmapped"

	if [ ${REG2} == 1 ]; then
		FLAGS="${FLAGS} --src_regional"
	fi

	if [ ${REG1} == 1 ]; then
		FLAGS="${FLAGS} --dst_regional"
	fi

	MAP_NAME="map_${GRID2_SHORTNAME}_TO_${GRID1_SHORTNAME}_${MAP_ABBREV}.${DATESTAMP}.nc"
	ESMF_RegridWeightGen --source ${GRID2} --destination ${GRID1} --method ${TYPE} --weight ${MAP_NAME} ${FLAGS} 
done
