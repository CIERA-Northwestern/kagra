#!/bin/bash
module load python
module load mpi/openmpi-1.8.3-intel2015.0

source /projects/b1011/non-lsc/lscsoft-user-env.sh
source /projects/b1011/ligo_project/lsc/lalinference_o2/etc/lscsoftrc

injxml=/projects/b1011/spinning_runs/IMRfreezinginj.xml

# NSBH: HLV vs HKLV
./full_sky_plot --underplot network \
    -p "HKLV=/projects/b1011/kagra/kagra_o2_lalinference/skymap/HKLV/*/post/*.fits.gz" -c "HKLV=error_region,hsv" \
    -p "HLV=/projects/b1011/kagra/kagra_o2_lalinference/skymap/HLV/*/post/*.fits.gz" -c "HLV=error_region,hsv" \
    --inj-xml $injxml --file-stem HLV_HKLV_comp

# NSBH: HKLV vs HIKLV
./full_sky_plot --underplot network \
    -p "HKLV=/projects/b1011/kagra/kagra_o2_lalinference/skymap/HLVK/*/post/*.fits.gz" -c "HKLV=error_region,hsv" \
    -p "HIKLV=/projects/b1011/kagra/kagra_o2_lalinference/skymap/HLVKI/*/post/*.fits.gz" -c "HIKLV=error_region,hsv" \
    --inj-xml $injxml --file-stem HKLV_HIKLV_comp
