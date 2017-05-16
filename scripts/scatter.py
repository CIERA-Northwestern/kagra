# Adapted from code developed by Carl Rodriguez

import glob
import os
from collections import defaultdict
from argparse import ArgumentParser
import operator

import numpy as np

from mpl_toolkits import basemap
from mpl_toolkits.basemap import Basemap
import matplotlib
matplotlib.use("agg")
from matplotlib import pyplot as plt
from matplotlib import cm as mpl_cm
from matplotlib import lines as mpl_lines
from matplotlib import colorbar
from matplotlib.ticker import MaxNLocator
import matplotlib.gridspec as gridspec

from skyutils import *

def get_event_num(fname):
    import re
    m = re.findall("([0-9]+)", fname)
    #print fname
    #print m
    return int(m[-1])

def load_injections(fname):
    from glue.ligolw import lsctables, utils, ligolw
    lsctables.use_in(ligolw.LIGOLWContentHandler)
    xmldoc = utils.load_filename(fname, contenthandler=ligolw.LIGOLWContentHandler)
    return lsctables.SimInspiralTable.get_table(xmldoc)

def parse_pathspecs(pspecs):
    out = {}
    for ps in pspecs:
        label, cmap = ps.split("=")
        try:
            cm = plt.get_cmap(cmap)
        except ValueError:
            cm = cmap
        out[label] = cmap
    return out

class Plottable(object):
    _scales = {
        "snr": {"vmin": 10, "vmax": 20},
        "error_region": {"vmin": 0, "vmax": 64},
        None: {"vmin": 0, "vmax": 64}
    }

    def __init__(self, colors, quant=None):
        self.colors = colors
        self.quant = quant
        self._init_cmaps()

    def _init_cmaps(self):
        try:
            self.colors = plt.get_cmap(self.colors)
            self._norm = \
                matplotlib.colors.Normalize(**Plottable._scales[self.quant])
        except ValueError:
            self._norm = None

    def __call__(self, val):
        if isinstance(self.colors, str):
            return self.colors
        return self.colors(self._norm(val))

    def __repr__(self):
        if isinstance(self.colors, str):
            return self.colors
        return self.colors.name + \
            " (%f, %f)" % (self._norm.vmin, self._norm.vmax)

def parse_colorspecs(pspecs):
    out = {}
    for ps in pspecs:
        label, cmap = ps.split("=")
        try:
            quant, cmap = cmap.split(",")
            out[label] = Plottable(cmap, quant)
        except ValueError:
            out[label] = Plottable(cmap)
    return out

_antenna_functions = ("network", "alignment", "dpf")
#_network = ("H1", "K1", "L1", "V1")
_network = "H1 K1 L1 V1"

argp = ArgumentParser()
argp.add_argument("-n", "--network", default=_network, help="Network of instruments to use for antenna based calculations. Default is %s" % ", ".join(_network))
argp.add_argument("-u", "--underplot", default=None, help="Underplot a function of the sky location as a continuous gradient. Valid choices are %s" % ", ".join(_antenna_functions))
argp.add_argument("-o", "--overplot", default=None, help="Overplot a function of the sky location as discrete contours. Valid choices are %s" % ", ".join(_antenna_functions))
argp.add_argument("-i", "--inj-xml", default=None, help="Path to the injection XML file.")
argp.add_argument("-p", "--pathspecs", action="append", help="Add file glob paths to be parsed in the following way: \"(name)=(globspec)\", e.g. \"HLV=/path/to/skymaps/*/post/*\"")
argp.add_argument("-c", "--colorspecs", action="append", help="Add color specs to be parsed in the following way: \"(name)=(colorspec)\", e.g. \"HLV=yellow\" or \"HKLV=snr,viridis\"")
argp.add_argument("-x", "--configuration", default=None, help="Choices are HLV, HKLV, HIKLV")
args = argp.parse_args()

network = args.network.split()
configuration =  args.configuration
#print "Will use network of %s" % ", ".join(network)

np.seterr(under='ignore')

plt.figure(0)
m = Basemap(projection='moll', lon_0=0.0, lat_0=0.0)

m.drawmapboundary()
m.drawcoastlines(linewidth=0.5)
m.drawparallels(np.arange(-90., 90., 45.), labels=[1, 0, 0, 0],
                labelstyle='+/-')
m.drawmeridians(np.arange(0., 360., 90.))  #,labels=[0,0,0,1],labelstyle='+/-')
net = "".join([n[0] for n in sorted(network)])
plt.title("Sky Location Credible Intervals, %s Configuration" % net,
          fontsize=20)
m.drawmapboundary()

# Arbitrary midnight UTC
gpstime = 1e9 - (1e9 % (3600 * 24))
ra_grid, dec_grid, net_pat, net_align, dpf = \
    net_antenna_pattern(gpstime, network) #, norm=True)
#print np.max(net_pat), np.max(net_align)

# Find Integrated network antenna pattern
filename = 'integrated_net_pat_%s' % configuration
integrated_net_pat = open(filename, 'w')
integrated_net_pat.write(str(integrate_net_pat(gpstime, network)))
integrated_net_pat.close()




#FIXME: Still getting warnings, not sure how this is supposed to work.
#net_pat, ra_grid = basemap.shiftgrid(180, net_pat.flatten(), ra_grid.flatten(), start=False)
#net_pat_tmp, ra_grid_tmp = basemap.shiftgrid(180., net_pat.flatten(), ra_grid.flatten())
#ra_grid_tmp = ra_grid_tmp.reshape(dec_grid.shape)
#net_pat_tmp = net_pat_tmp.reshape(dec_grid.shape)

# Probably because we're viewing the sphere from "outside in".
ra_grid -= 180
dec_grid *= -1

"""
# TODO: Make into a function and be consistent
ra_grid = ra_grid * 180 / np.pi
ra_grid = np.where(ra_grid < -180, ra_grid + 360, ra_grid)
ra_grid = np.where(ra_grid > 180, ra_grid - 360, ra_grid)
# ANNOYING!
while ra_grid[0,0]/ra_grid[0,-1] > 0:
    ra_grid = np.roll(ra_grid, -1, axis=1)
    net_pat = np.roll(net_pat, -1, axis=1)
    net_align = np.roll(net_align, -1, axis=1)
    dpf = np.roll(dpf, -1, axis=1)
ra_grid *= -1
"""

ra_grid, dec_grid = m(ra_grid, dec_grid)

#
# Underplot various network quantities
#
if args.underplot == "network":
    m.contourf(ra_grid, dec_grid, net_pat, 1000, cmap=matplotlib.cm.Greys_r)
    #plt.savefig("net_pat.png")

elif args.underplot == "dpf":
    m.contourf(ra_grid, dec_grid, dpf, 1000, cmap=matplotlib.cm.hsv,
                vmin=-np.pi/4, vmax=np.pi/4)
    #plt.savefig("net_dpf.png")

elif args.underplot == "align":
    m.contourf(ra_grid, dec_grid, net_align, 1000, cmap=matplotlib.cm.Greys_r, vmin=0, vmax=1)
    #plt.savefig("net_align.png")

#
# Overplot contours of the network antenna pattern
#
if args.overplot == "network":
    m.contour(ra_grid, dec_grid, net_pat, 10, cmap=matplotlib.cm.cool)
    #plt.savefig("net_pat.png")

elif args.overplot == "dpf":
    m.contour(ra_grid, dec_grid, dpf, 10, cmap=matplotlib.cm.cool,
                vmin=-np.pi/4, vmax=np.pi/4)
    #plt.savefig("net_dpf.png")

elif args.overplot == "align":
    m.contour(ra_grid, dec_grid, net_align, 10, cmap=matplotlib.cm.cool,
                vmin=0, vmax=1)
    #plt.savefig("net_align.png")

#
# Load injections
#
if args.inj_xml is not None:
    inj = load_injections(args.inj_xml)
    #print "Loaded %d injections" % len(inj)

cspecs = parse_colorspecs(args.colorspecs)
#for label, cmap in cspecs.iteritems():
    #print label + " => " + str(cmap)

config_information = defaultdict(list)
outliers = []
pspecs = parse_pathspecs(args.pathspecs)
for label, globpat in pspecs.iteritems():

    #print label + " => " + globpat

    files = glob.glob(globpat)
    #print "Globbed %d files for pattern %s" % (len(files), globpat)

    plt.figure(0)
    for filename in files: #[:3]:  # Change this to run faster
        enum = get_event_num(filename)
        #print "Processing event %d" % enum

        #
        # Generate vital statistics
        #
        gmst = inj[enum].geocent_end_time
        gmst = np.mod(gmst/3600., 24.)

        sky_data, smap = create_table(filename)
        ns = healpy.npix2nside(len(smap))
        pix_size = healpy.nside2pixarea(ns, degrees=True)

        prb68 = np.searchsorted(sky_data["cumul"], 0.68)
        prb90 = np.searchsorted(sky_data["cumul"], 0.90)
        prb95 = np.searchsorted(sky_data["cumul"], 0.95)
        prb99 = np.searchsorted(sky_data["cumul"], 0.99)

        if configuration == "HLV":
            snr = "/projects/b1011/spinning_runs/freezingparams_20160402_IMR/" + str(enum)  + "/none/snr.txt"
        elif configuration == "HKLV":
            snr = "/projects/b1011/kagra/kagra_o2_lalinference/" + str(enum)  + "/snr.txt"
        elif configuration == "HIKLV":
            snr = "/projects/b1011/kagra/HLVKI/" + str(enum) + "/snr.txt"


        #print snr
        try:
            with open(snr, "r") as snrf:
                snrs = dict([map(str.strip, l.split(":")) for l in snrf.readlines()])
            for k in snrs:
                snrs[k] = float(snrs[k])
            #print snrs
            # Make a dictionary without the Network value
            snrs_reduced = dict(snrs)
            del snrs_reduced['Network']
            print snrs_reduced

            # Remove the detector with the highest SNR
            detector_high_snr = sorted(snrs_reduced.iteritems(), key=operator.itemgetter(1))[-1][0]
            del snrs_reduced[detector_high_snr]
            print snrs_reduced

            # Find the detector with the second-highest SNR
            detector_second_snr = sorted(snrs_reduced.iteritems(), key=operator.itemgetter(1))[-1][0]
            # FIXME: change the threshold to be 5.5 and also check that the new dict sorting is working
            # Check value of second-highest SNR
            if snrs[detector_second_snr] < 5.0:
                #outliers.append((configuration, enum, snrs[detector_second_snr]))
                outliers.append(enum)
                print snrs
            # Add the special case, 901, to the outlier list
            if enum == 901:
                outliers.append(enum)

        except IOError:
            snrs_new = {}
            snrs_new["Network"] = 1.0
            snrs = snrs_new

        #print snrs["Network"]

        if False:
            cmap = cspecs[label]
            if cmap.quant == "error_regior":
                linecolor = cmap(prb90 * pix_size)
            elif cmap.quant == "snr":
                linecolor = cmap(snrs["Network"])

            cache_fname = filename.replace("fits.gz", "intrp.npz")
            if os.path.exists(cache_fname):
                #print "Loading cached interpolated map: %s" % cache_fname
                ra_int, dec_int, prob_int = dict(np.load(cache_fname))["arr_0"]
            else:
                ra_int, dec_int, prob_int = interpolate_healpix_map(sky_data["ra"],
                    sky_data["dec"], smap, npts=500)
                np.savez(cache_fname, [ra_int, dec_int, prob_int])

            # Better idea, figure out by how much this shift the ra_int and shift
            # the z array by the opposite of that amount
            ra_int = (-gmst + ra_int) * np.pi / 12
            # Reverse... and undo the radians...
            #ra_int = -ra_int * 180 / np.pi
            ra_int = ra_int * 180 / np.pi
            ra_int = np.where(ra_int < -180, ra_int + 360, ra_int)
            ra_int = np.where(ra_int > 180, ra_int - 360, ra_int)

            # ANNOYING!
            while ra_int[0,0]/ra_int[0,-1] > 0:
                ra_int = np.roll(ra_int, -1, axis=1)
                prob_int = np.roll(prob_int, -1, axis=1)
            ra_int *= -1

            ra_int, dec_int = m(ra_int, dec_int)

            m.contour(ra_int, dec_int, prob_int, [sky_data["prob"][prb90]], colors=(linecolor,), linewidths=0.5)

        # Use gpstime of this injection to find the network antenna pattern at that time at this ra and dec
        antenna_pattern = net_antenna_pattern_point(gmst, network, inj[enum].longitude, inj[enum].latitude)[0]
        alignment = net_antenna_pattern_point(gmst, network, inj[enum].longitude, inj[enum].latitude)[1]
        config_information[label].append([prb90 * pix_size, snrs["Network"], antenna_pattern, enum, alignment])

        #print gmst, network, inj[enum].longitude, inj[enum].latitude
            # Debuggin
            #m.scatter(ra_int.flatten()[-200000:], dec_int.flatten()[-200000:], c=prob_int.flatten()[-200000:], marker='.', edgecolor='none')

        # FIXME: temporary
        plt.savefig("figures/test_map.png")

#filename = 'outliers_%s' % configuration
#outlier_file = open(filename, 'w')
#for run in outliers:
    #outlier_file.write('Event number: ' + str(run[1]) + ', Second-Highest SNR: ' + str(run[2]))
    #outlier_file.write('\n')
#    outlier_file.write(str(run))
#    outlier_file.write('\n')

#outlier_file.write(outliers)
#outlier_file.close()
#print outliers

plt.figure(1)
gs = gridspec.GridSpec(3, 3)
gs.update(hspace=0.0, wspace=0.0)
ax1 = plt.subplot(gs[1:, :-1])
ax2 = plt.subplot(gs[0, :-1], sharex=ax1)
ax3 = plt.subplot(gs[1:, 2], sharey=ax1)

filename = 'err_snr_antenna_%s' % configuration
plot_data = open(filename, 'w')


for label, config in config_information.iteritems():
    config = np.asarray(config)
    err_reg = config[:, 0]
    snr = config[:, 1]
    antenna_pattern = config[:, 2]

    for value in config:
        plot_data.write(str(value[0]) + ' ' + str(value[1]) + ' ' + str(value[2]) + ' ' + str(int(value[3])) + ' ' + str(value[4]))
        plot_data.write('\n')

    # Scatter plot of SNR vs. Error Regions
    scatter = ax1.scatter(err_reg, antenna_pattern, c=snr, edgecolors='None', label=label, cmap='viridis')
    #ax1.tick_params(top='off', bottom='on', left='on', right='off')
    ax1.set_xlabel('Error Region (squared degrees)')
    ax1.set_ylabel('Network Antenna Pattern')
    ax1.set_xlim(1e-1, 1e5)
    ax1.set_ylim(0, 4)
    ax1.set_xscale('log')

    # Histogram of Error Regions
    ax2.hist(err_reg, color='gray', histtype='stepfilled', bins=np.logspace(-1, 5, 20))
    ax2.tick_params(axis='x', top='on', bottom='on', labelbottom='off', labeltop='off')
    ax2.set_ylabel('Count')
    nbins = len(ax1.get_xticklabels())
    ax2.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='lower', integer=True))

    # Histogram of SNRs
    ax3.hist(antenna_pattern, color='gray', histtype='stepfilled', orientation='horizontal')
    ax3.tick_params('y', left='on', right='on', labelleft='off', labelright='off')
    ax3.set_xlabel('Count')
    #nbins = len(ax3.get_xticklabels())
    ax3.xaxis.set_major_locator(MaxNLocator(integer=True))
    xticks = ax3.xaxis.get_major_ticks()
    xticks[0].label1.set_visible(False)


    #ax3.set_yticklabels([])

    #ax1.xlabel('Error Region (squared degrees)')
    #ax1.ylabel('Network SNR')
    plt.suptitle('SNR, Error Regions, and Network Antenna Pattern for %s' % configuration)
    #plt.subplots_adjust(right=0.8)
    #cbar_ax = plt.add_axes([0.85, 0.15, 0.05, 0.7])
    cbar = plt.colorbar(scatter, label='Network SNR')
    cbar.ax.tick_params(labelsize=10)


    #ax1.grid()
    #plt.step(err_reg[0], yaxis, label=label)

plt.savefig("figures/snr_vs_err_%s.png" % configuration)
plot_data.close()
