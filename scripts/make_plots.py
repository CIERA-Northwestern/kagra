import glob
import os
from collections import defaultdict
from argparse import ArgumentParser

import numpy as np

from mpl_toolkits import basemap
from mpl_toolkits.basemap import Basemap
import matplotlib
matplotlib.use("agg")
from matplotlib import pyplot as plt
from matplotlib import cm as mpl_cm
from matplotlib import lines as mpl_lines
from matplotlib import colorbar
from matplotlib import colors
from matplotlib.ticker import MaxNLocator
import matplotlib.gridspec as gridspec

from skyutils import *

# Determine which plots to make
argp = ArgumentParser()
argp.add_argument("-n", "--network", default=None, help="Choices are all, HLV, 'HLV_to_HKLV', and 'HKLV_to_HIKLV'")
argp.add_argument("-y", "--yaxis", default='net_pat', help="Choices are 'net_pat', 'align', and 'angle_err'")
argp.add_argument("-i", "--hist", default='no_net_hist', help="Choices are 'no_net_hist' and 'net_hist'")
args = argp.parse_args()

histogram_type = args.hist
if histogram_type == 'no_net_hist':
    net_hist = False
elif histogram_type == 'net_hist':
    net_hist = True

network_to_plot = args.network

runs_kagra = np.loadtxt('err_snr_antenna_HKLV')
runs_nokagra = np.loadtxt('err_snr_antenna_HLV')
runs_india = np.loadtxt('err_snr_antenna_HIKLV')

outlier_list = np.loadtxt('outliers')

if network_to_plot == 'all':
    runs_to_plot = [runs_nokagra, runs_kagra, runs_india]
    names = ['HLV', 'HKLV', 'HIKLV']
    scatter_config = 'HIKLV'
elif network_to_plot == 'HLV_to_HKLV':
    runs_to_plot = [runs_nokagra, runs_kagra]
    names = ['HLV', 'HKLV']
    scatter_config = 'HKLV'
elif network_to_plot == 'HKLV_to_HIKLV':
    runs_to_plot = [runs_kagra, runs_india]
    names = ['HKLV', 'HIKLV']
    scatter_config = 'HIKLV'
elif network_to_plot == 'HLV':
    runs_to_plot = [runs_nokagra]
    names = ['HLV']
    scatter_config = 'HLV'



# Set up plots
plt.figure(1)
#gs = gridspec.GridSpec(3, 3)
#gs.update(hspace=0.0, wspace=0.0)
if net_hist:
    gs = gridspec.GridSpec(3, 3)
    gs.update(hspace=0.0, wspace=0.0)
    ax1 = plt.subplot(gs[1:, :-1])
    ax2 = plt.subplot(gs[0, :-1], sharex=ax1)
    #if net_hist:
    ax3 = plt.subplot(gs[1:, 2], sharey=ax1)
else:
    gs = gridspec.GridSpec(30, 30)
    gs.update(hspace=0.0, wspace=0.0)
    ax1 = plt.subplot(gs[10:, :-2])
    ax2 = plt.subplot(gs[0:10, :-2], sharex=ax1)
    ax3 = plt.subplot(gs[10:, 29]) #, sharey=ax1)

#print outlier_list
#print outlier_list[0]

scatter_shape = ['o', '+', 's']
colors_options = ['gray', 'cyan', 'white']
colors_stats = ['black', 'blue', 'red']

# Print median and 90% interval of err regions to file
if network_to_plot == 'all':
    interval_file = open('median_and_90conf', 'w')


for runs in runs_to_plot:

    #print runs[0, :]

    err_reg = runs[:, 0]
    snr = runs[:, 1]
    antenna_pattern = runs[:, 2]
    run_numbers = runs[:, 3]
    alignment = runs[:, 4]
    angle_err = runs[:, 5]
    hist_color = colors_options.pop(0)
    statistics_color = colors_stats.pop(0)
    label = names.pop(0)

    yvariable = args.yaxis

    if yvariable == 'net_pat':
        yaxis = antenna_pattern
    elif yvariable == 'align':
        yaxis = alignment
    elif yvariable == 'angle_err':
        yaxis = angle_err

    # Find points that aren't outliers for histograms
    err_reg_noout = []
    yaxis_noout = []

    # FIXME: list comprehension
    for index, enum in enumerate(run_numbers):
        if enum not in outlier_list:
            err_reg_noout.append(err_reg[index])
            yaxis_noout.append(yaxis[index])


    #Scatter plot of SNR vs. Error Regions
    if label == scatter_config:

        norm = colors.Normalize(vmin=min(snr), vmax=max(snr))
        smap = mpl_cm.ScalarMappable(norm=norm, cmap=plt.get_cmap("viridis"))

        nonoutlier_marker = scatter_shape.pop(0)
        outlier_marker = scatter_shape.pop(0)
        first_plotted_point = True
        for index, enum in enumerate(run_numbers):
            if enum not in outlier_list:
                marker = nonoutlier_marker
                if first_plotted_point:
                    scatter = ax1.scatter(err_reg[index], yaxis[index], marker=marker, c=smap.to_rgba(snr[index]), edgecolors='None', cmap='viridis', label=scatter_config)
                    first_plotted_point = False
                else:
                    scatter = ax1.scatter(err_reg[index], yaxis[index], marker=marker, c=smap.to_rgba(snr[index]), edgecolors='None', cmap='viridis')



        if yvariable == 'net_pat':
            ax1.set_ylim(0.0, 2)
            ax1.set_ylabel('Network Antenna Pattern')

        elif yvariable == 'align':
            ax1.set_ylim(0.0, 1.0)
            ax1.set_ylabel('Network Alignment Factor')

        elif yvariable == 'angle_err':
            ax1.set_ylim(0.0, None)
            ax1.set_ylabel('Solid Angle Error')


        ax1.set_xlabel('Error Region (squared degrees)')
        ax1.set_xlim(1e-1, 1e5)
        ax1.set_xscale('log')

    # Histogram of Error Regions
    ax2.hist(err_reg_noout, color=hist_color, histtype='stepfilled', bins=np.logspace(-1, 5, 20), alpha=0.5, label=label)
    ax2.tick_params(axis='x', top='on', bottom='on', labelbottom='off', labeltop='off')
    #ax2.set_ylabel('Count')
    nbins = len(ax1.get_xticklabels())
    ax2.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='lower', integer=True))
    ax2.get_yaxis().set_ticks([])


    # Find the median and 90% interval of Error Region histograms
    err_med = np.percentile(err_reg_noout, 50)
    err_90 = np.percentile(err_reg_noout, 90)
    #err_5 = np.percentile(err_reg_noout, 5)
    ax2.axvline(x=err_med, color=statistics_color, linewidth=2, label=(label + ' Median'))
    ax2.axvline(x=err_90, color=statistics_color, linewidth=2, ls='-.', label=(label + ' 90% Conf.'))
    #ax2.axvline(x=err_, color=statistics_color, linewidth=2, ls='-.')

    # Write the 50% and 90% regions to file
    if network_to_plot == 'all':
        interval_file.write('Network: %s, Median: %.2f, 90 Percent Confidence Interval: %.2f' % (
            label, err_med, err_90))
        interval_file.write('\n')

    # Histogram of Antenna Patterns
    if yvariable == 'net_pat':
        bins = np.linspace(0, 2, 20)
    elif yvariable == 'align':
        bins = np.linspace(0, 1, 20)
    elif yvariable == 'angle_err':
        bins = np.linspace(0, max(yaxis_noout), 20)

    if net_hist:
        ax3.hist(yaxis_noout, color=hist_color, histtype='stepfilled', bins=bins, orientation='horizontal', alpha=0.5)
        ax3.tick_params('y', left='on', right='on', labelleft='off', labelright='off')
        #ax3.set_xlabel('Count')
        ax3.xaxis.set_major_locator(MaxNLocator(integer=True))
        xticks = ax3.xaxis.get_major_ticks()
        xticks[0].label1.set_visible(False)
        ax3.get_xaxis().set_ticks([])


# Close interval file
if network_to_plot == 'all':
    interval_file.close()

# Overall plot features
smap.set_array([0, 1])
if net_hist:
    cbar = plt.colorbar(smap, label='Network SNR')
else:
    cbar = plt.colorbar(smap, label='Network SNR', cax=ax3) #ax=ax1, orientation='horizontal')
cbar.ax.tick_params(labelsize=7)
ax1.legend(loc='upper right', fontsize=7)
ax2.legend(loc='upper right', fontsize=7)
#ax3.legend(loc='lower right', fontsize=7)

plt.savefig("figures/%s_%s.png" % (network_to_plot, yvariable))
