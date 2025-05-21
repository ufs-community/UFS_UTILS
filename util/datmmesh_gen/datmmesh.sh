#!/bin/bash
set -eux

APRUN=${APRUN:-"srun"}

FDIMS=${NX}x${NY}
FDST=${OUTPUT_DIR}/datm.${FDIMS}.SCRIP.nc
if [ $N2S == .true. ]; then
    ncremap -g ${FDST} -G ttl='DATM grid '${FDIMS}#latlon=${NY},${NX}#lon_typ=grn_ctr#lat_typ=gss#lat_drc=n2s
else
    ncremap -g ${FDST} -G ttl='DATM grid '${FDIMS}#latlon=${NY},${NX}#lon_typ=grn_ctr#lat_typ=gss
fi

FSRC=${OUTPUT_DIR}/datm.${FDIMS}.SCRIP.nc
FDST=${OUTPUT_DIR}/mesh.datm.${FDIMS}.nc
$APRUN -n 1 ESMF_Scrip2Unstruct ${FSRC} ${FDST} 0
