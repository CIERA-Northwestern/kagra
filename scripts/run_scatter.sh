#!/bin/bash
module load python
module load mpi/openmpi-1.8.3-intel2015.0

source /projects/b1011/non-lsc/lscsoft-user-env.sh
source /projects/b1011/ligo_project/lsc/lalinference_o2/etc/lscsoftrc

# FIXME: several of these labels are for "fixed" runs without actually pointing to fixed runs

# Fixed injection HKLV -- change this
#python scatter.py --underplot network -p "HLVK_fixed=/projects/b1011/kagra/HLVK/skymap/*/*.fits.gz" -c "HLVK_fixed=snr,hsv" --inj-xml /projects/b1011/spinning_runs/IMRfreezinginj.xml -x "HKLV"

# Fixed injection HLV (no Kagra)
python scatter.py --underplot network -p "HLV_fixed=/projects/b1011/kagra/HLV/*/none/post/*.fits.gz" -c "HLV_fixed=snr,hsv" --inj-xml /projects/b1011/spinning_runs/IMRfreezinginj.xml -x "HLV" -n "H1 L1 V1"

# Fixed injection HILKV
#python scatter.py --underplot network -p "HLVKI_fixed=/projects/b1011/kagra/HLVKI/skymap/*/*.fits.gz" -c "HLVKI_fixed=snr,hsv" --inj-xml /projects/b1011/spinning_runs/IMRfreezinginj.xml -x "HIKLV" -n "H1 L1 V1 K1 I1"







