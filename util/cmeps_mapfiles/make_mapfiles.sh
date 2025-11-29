#!/bin/bash
set -eux

APRUN=${APRUN:-"srun -n 24"}

# Should be provided by the caller (e.g. `run.ursa.sh`).
ATMRES=${ATMRES:-}
OCNRES=${OCNRES:-}
WAVRES=${WAVRES:-}
FIX_DIR=${FIX_DIR:-}
orog_ver=${orog_ver:-}
ice_ver=${ice_ver:-}
wav_ver=${wav_ver:-}

fv3dir="${FIX_DIR}/orog/${orog_ver}"
datmdir=/scratch4/NAGAPE/epic/role-epic/UFS-WM_RT/NEMSfv3gfs/input-data-20251015/DATM_CDEPS
#icedir="${FIX_DIR}/cice/${ice_ver}"
wavdir="${FIX_DIR}/wav/${wav_ver}"

icedir=/scratch4/NCEPDEV/stmp/Denise.Worthen/CPLD_GRIDGEN/BASELINE

# Set ATM mesh based on ATMRES
if [[ $ATMRES == C* ]]; then
    # FV3 cube-sphere grid
    fmosaic="${fv3dir}/${ATMRES}/${ATMRES}_mosaic.nc"
    ftilepath="${fv3dir}/${ATMRES}"
    fatmmesh=""
else
    # DATM unstructured mesh
    fatmmesh="${datmdir}/mesh.datm.${ATMRES}.nc"
    fmosaic=""
    ftilepath=""
fi

# Set ocean mesh if OCNRES is provided
if [ -n "${OCNRES:-}" ]; then
    focnmesh=$icedir/${OCNRES}/'mesh.mx'${OCNRES}'.nc'
fi

# Set wave mesh if WAVRES is provided
if [ -n "${WAVRES:-}" ]; then
    fwavmesh="${wavdir}/mesh.${WAVRES}.nc"
fi

# Set srcopt based on ATM mesh type
if [ -n "${fmosaic}" ] && [ -n "${ftilepath}" ]; then
    srcopt="-s ${fmosaic} --tilefile_path ${ftilepath}"
elif [ -n "${fatmmesh}" ]; then
    srcopt="-s ${fatmmesh}"
else
    echo "Error: no ATM grid specified (set fmosaic+ftilepath or fatmmesh)" >&2
    exit 1
fi

defaultopts=' --src_loc center --dst_loc center --weight_only --no_log'
#defaultopts=' --src_loc center --dst_loc center --checkFlag '

for exp in a2o_bilin a2o_consf a2o_patch a2w_bilin w2o o2w; do
    # Skip ocean-related mappings if OCNRES not provided
    if [ -z "${OCNRES:-}" ] && [[ ${exp:0:3} == *o* ]]; then
        continue
    fi

    # Skip wave-related mappings if WAVRES not provided
    if [ -z "${WAVRES:-}" ] && [[ ${exp:0:3} == *w* ]]; then
        continue
    fi

    case $exp in
        w2o)
            mapindex=bilnr_nstod
            ftag=${OUTPUT_DIR}/'map.'${WAVRES}'.to.mx'${OCNRES}'.'$mapindex'.nc'
            mapping='-m bilinear -p none --extrap_method neareststod '
            opts='-s '${fwavmesh}' -d '${focnmesh}' -w '${ftag}'  '${mapping}
            ;;
        o2w)
            mapindex=bilnr_nstod
            ftag=${OUTPUT_DIR}/'map.mx'${OCNRES}'.to.'${WAVRES}'.'$mapindex'.nc'
            mapping='-m bilinear -p none --extrap_method neareststod '
            opts='-s '${focnmesh}' -d '${fwavmesh}' -w '${ftag}'  '${mapping}
            ;;
        a2o_bilin)
            mapindex=bilnr
            ftag=${OUTPUT_DIR}/'map.'${ATMRES}'.to.mx'${OCNRES}'.'$mapindex'.nc'
            mapping='-m bilinear -p all '
            opts="${srcopt} -d ${focnmesh} -w ${ftag} ${mapping}"
            ;;
        a2o_consf)
            mapindex=consf
            ftag=${OUTPUT_DIR}/'map.'${ATMRES}'.to.mx'${OCNRES}'.'$mapindex'.nc'
            mapping='-m conserve --norm_type fracarea '
            opts="${srcopt} -d ${focnmesh} -w ${ftag} ${mapping}"
            ;;
        a2o_patch)
            mapindex='patch'
            ftag=${OUTPUT_DIR}/'map.'${ATMRES}'.to.mx'${OCNRES}'.'$mapindex'.nc'
            mapping='-m patch -p all '
	    opts="${srcopt} -d ${focnmesh} -w ${ftag} ${mapping}"
            ;;
        a2w_bilin)
            mapindex=bilnr
            ftag=${OUTPUT_DIR}/'map.'${ATMRES}'.to.'${WAVRES}'.'$mapindex'.nc'
            mapping='-m bilinear -p none '
	    opts="${srcopt} -d ${fwavmesh} -w ${ftag} ${mapping}"
	    ;;
    esac

    ${APRUN} ESMF_RegridWeightGen ${opts} ${defaultopts}
done
