#!/bin/bash
module load python
module load mpi/openmpi-1.8.3-intel2015.0

source /projects/b1011/non-lsc/lscsoft-user-env.sh
source /projects/b1011/ligo_project/lsc/lalinference_o2/etc/lscsoftrc

# Fixed injection HLKV
python scatter.py --underplot network -p "HLVK_fixed=/projects/b1011/kagra/kagra_o2_lalinference/fixed_inj/*/post/*.fits.gz" -c "HLVK_fixed=snr,hsv" --inj-xml /projects/b1011/kagra/kagra_o2_lalinference/fixed_inj/fixed_inj.xml.gz

# Fixed injection HLV (no Kagra)
python scatter.py --underplot network -p "HLV_fixed=/projects/b1011/kagra/kagra_o2_lalinference/fixed_inj/HLV/fixed_distance/*/post/*.fits.gz" -c "HLV_fixed=error_region,hsv" --inj-xml /projects/b1011/kagra/kagra_o2_lalinference/fixed_inj/fixed_inj_hklv.xml.gz

