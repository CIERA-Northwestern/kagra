import numpy as np
from argparse import ArgumentParser

argp = ArgumentParser()
argp.add_argument("-n", "--network", default=None, help="Choices are HLV, HKLV, HIKLV")
args = argp.parse_args()
network = args.network
runs = np.loadtxt('err_snr_antenna_%s' % network)


# Remove all outliers with less than SNR 20
runs_reduced = filter(lambda x: x[1] >= 20, runs)

# Find largest error region out of the remaining values
largest_err = runs_reduced[np.argmax(runs_reduced, axis=0)[0]]
print 'Event over 20 SNR with largest error region: Event #%i with error region %.2f sq. deg. and SNR %.2f' % (
    largest_err[3], largest_err[0], largest_err[1])




