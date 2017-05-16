#!/bin/bash
module load python
module load mpi/openmpi-1.8.3-intel2015.0

source /projects/b1011/non-lsc/lscsoft-user-env.sh
source /projects/b1011/ligo_project/lsc/lalinference_o2/etc/lscsoftrc

injxml=/projects/b1011/spinning_runs/IMRfreezinginj.xml

# NSBH: HLV vs HKLV
/bin/true ./full_sky_plot --underplot network \
    -p "HKLV=/projects/b1011/kagra/HLVK/skymap/*/*.fits.gz" -c "HKLV=error_region,hsv" \
    -p "HLV=/projects/b1011/kagra/kagra_o2_lalinference/skymap/HLV/*/none/post/*.fits.gz" -c "HLV=error_region,hsv" \
    --inj-xml $injxml --file-stem HLV_HKLV_comp --black-list outliers.txt

./full_sky_plot --underplot network \
    -p "HKLV=/projects/b1011/kagra/HLVK/skymap/*/*.fits.gz" -c "HKLV=error_region,green" \
    -p "HLV=/projects/b1011/kagra/kagra_o2_lalinference/skymap/HLV/*/none/post/*.fits.gz" -c "HLV=error_region,red" \
    --inj-xml $injxml --file-stem HLV_HKLV_comp --black-list outliers.txt

# NSBH: HKLV vs HIKLV
./full_sky_plot --underplot network --network "H1 I1 K1 L1 V1" \
    -p "HKLV=/projects/b1011/kagra/kagra_o2_lalinference/skymap/HLVK/*/post/*.fits.gz" -c "HKLV=error_region,hsv" \
    -p "HIKLV=/projects/b1011/kagra/kagra_o2_lalinference/skymap/HLVKI/*/*.fits.gz" -c "HIKLV=error_region,hsv" \
    --inj-xml $injxml --file-stem HKLV_HIKLV_comp --black-list outliers.txt

./full_sky_plot --underplot network --network "H1 I1 K1 L1 V1" \
    -p "HKLV=/projects/b1011/kagra/kagra_o2_lalinference/skymap/HLVK/*/post/*.fits.gz" -c "HKLV=error_region,green" \
    -p "HIKLV=/projects/b1011/kagra/kagra_o2_lalinference/skymap/HLVKI/*/*.fits.gz" -c "HIKLV=error_region,cyan" \
    --inj-xml $injxml --file-stem HKLV_HIKLV_comp --black-list outliers.txt
