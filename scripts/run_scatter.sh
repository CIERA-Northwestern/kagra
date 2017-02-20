#!/bin/bash
module load python
module load mpi/openmpi-1.8.3-intel2015.0

source /projects/b1011/non-lsc/lscsoft-user-env.sh
source /projects/b1011/ligo_project/lsc/lalinference_o2/etc/lscsoftrc

# FIXME: several of these labels are for "fixed" runs without actually pointing to fixed runs

# Fixed injection HLKV --- FIXME: you are not using the correct fitz.gz file. change where we look for SNR
python scatter.py --underplot network -p "HLVK_fixed=/projects/b1011/kagra/kagra_o2_lalinference/skymap/HLVK/*/post/*.fits.gz" -c "HLVK_fixed=snr,hsv" --inj-xml /projects/b1011/spinning_runs/IMRfreezinginj.xml -x "HLKV"

# Fixed injection HLV (no Kagra)
python scatter.py --underplot network -p "HLV_fixed=/projects/b1011/kagra/kagra_o2_lalinference/skymap/HLV/*/none/post/*.fits.gz" -c "HLV_fixed=snr,hsv" --inj-xml /projects/b1011/spinning_runs/IMRfreezinginj.xml -x "HLV"

# Normalize error region: divide by sqrt(2) but not by the detectors







