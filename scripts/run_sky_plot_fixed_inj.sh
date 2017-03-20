#!/bin/bash
module load python
module load mpi/openmpi-1.8.3-intel2015.0

source /projects/b1011/non-lsc/lscsoft-user-env.sh
source /projects/b1011/ligo_project/lsc/lalinference_o2/etc/lscsoftrc

# Fixed injection run1
#./full_sky_plot --underplot network -p "HLVK_fixed=/projects/b1011/kagra/kagra_o2_lalinference/fixed_inj/run1/*/post/*.fits.gz" -c "HLVK_fixed=snr,hsv" --inj-xml /projects/b1011/kagra/kagra_o2_lalinference/fixed_inj/fixed_inj.xml.gz --file-stem HKVL_snr
#./full_sky_plot --underplot network -p "HLVK_fixed=/projects/b1011/kagra/kagra_o2_lalinference/fixed_inj/run1/*/post/*.fits.gz" -c "HLVK_fixed=error_region,hsv" --inj-xml /projects/b1011/kagra/kagra_o2_lalinference/fixed_inj/fixed_inj.xml.gz --file-stem HKVL_err

# Fixed injection HLV vs HKLV, fixed distance
./full_sky_plot --underplot network \
    -p "HKLV_fixed=/projects/b1011/kagra/kagra_o2_lalinference/fixed_inj/HKLV/fixed_distance/*/post/*.fits.gz" -c "HKLV_fixed=error_region,hsv" \
    -p "HLV_fixed=/projects/b1011/kagra/kagra_o2_lalinference/fixed_inj/HLV/fixed_distance/*/post/*.fits.gz" -c "HLV_fixed=error_region,hsv" \
    --inj-xml /projects/b1011/kagra/kagra_o2_lalinference/fixed_inj/fixed_inj_hklv.xml.gz --file-stem fixed_dist_comp

# Unfixed injection HLV vs HKLV, fixed distance
./full_sky_plot --underplot network \
    -p "HKLV_unfixed=/projects/b1011/kagra/kagra_o2_lalinference/fixed_inj/HKLV/unfixed_distance/*/post/*.fits.gz" -c "HKLV_unfixed=error_region,hsv" \
    -p "HLV_unfixed=/projects/b1011/kagra/kagra_o2_lalinference/fixed_inj/HLV/unfixed_distance/*/post/*.fits.gz" -c "HLV_unfixed=error_region,hsv" \
    --inj-xml /projects/b1011/kagra/kagra_o2_lalinference/fixed_inj/fixed_inj_hklv.xml.gz --file-stem unfixed_dist_comp

# Fixed vs unfixed injection HLV
./full_sky_plot --underplot network \
    -p "HLV_fixed=/projects/b1011/kagra/kagra_o2_lalinference/fixed_inj/HLV/fixed_distance/*/post/*.fits.gz" -c "HLV_fixed=error_region,hsv" \
    -p "HLV_unfixed=/projects/b1011/kagra/kagra_o2_lalinference/fixed_inj/HLV/unfixed_distance/*/post/*.fits.gz" -c "HLV_unfixed=error_region,hsv" \
    --inj-xml /projects/b1011/kagra/kagra_o2_lalinference/fixed_inj/fixed_inj_hlv.xml.gz --file-stem hlv_comp

# Fixed vs unfixed injection HKLV
./full_sky_plot --underplot network \
    -p "HKLV_fixed=/projects/b1011/kagra/kagra_o2_lalinference/fixed_inj/HKLV/fixed_distance/*/post/*.fits.gz" -c "HKLV_fixed=error_region,hsv" \
    -p "HKLV_unfixed=/projects/b1011/kagra/kagra_o2_lalinference/fixed_inj/HKLV/unfixed_distance/*/post/*.fits.gz" -c "HKLV_unfixed=error_region,hsv" \
    --inj-xml /projects/b1011/kagra/kagra_o2_lalinference/fixed_inj/fixed_inj_hklv.xml.gz --file-stem hklv_comp
