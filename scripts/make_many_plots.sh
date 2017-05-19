#!/bin/bash
module load python
module load mpi/openmpi-1.8.3-intel2015.0

source /projects/b1011/non-lsc/lscsoft-user-env.sh
source /projects/b1011/ligo_project/lsc/lalinference_o2/etc/lscsoftrc


# Make plots for HLV to HKLV
python make_plots.py --network 'HLV_to_HKLV' -y 'angle_err'

exit

# Make plots with all configurations
python make_plots.py --network 'all' -y 'align'

# Make plots for HLV to HKLV
python make_plots.py --network 'HLV_to_HKLV' -y 'align'

# Make plots for HKLV to HIKLV
python make_plots.py --network 'HKLV_to_HIKLV' -y 'align'

# Just HLV
python make_plots.py --network 'HLV' -y 'align'


# Plots for network antenna pattern
# Make plots with all configurations
python make_plots.py --network 'all'

# Make plots for HLV to HKLV
python make_plots.py --network 'HLV_to_HKLV'

# Make plots for HKLV to HIKLV
python make_plots.py --network 'HKLV_to_HIKLV'

# Just HLV
python make_plots.py --network 'HLV'





