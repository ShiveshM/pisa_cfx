#!/usr/bin/env python
"""
From a MapSet json file, create an efficiency plot for all flavours
"""

import os
import argparse
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

import numpy as np
from uncertainties import unumpy as unp
import matplotlib as mpl
from matplotlib import pyplot as plt
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
mpl.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
mpl.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'

from pisa.core.map import MapSet
from pisa.utils import fileio
from pisa.utils.log import logging, set_verbosity

class FullPaths(argparse.Action):
    """
    Append user- and relative-paths
    """
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))

def parse_args():
    """Get command line arguments"""
    parser = ArgumentParser(
        description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '-i', '--infile', type=fileio.is_valid_file, action=FullPaths,
        metavar='FILE',
        help='location/filename of the json file containing the MapSet'
    )
    parser.add_argument(
        '-o', '--outdir', type=fileio.is_dir, action=FullPaths, metavar='DIR',
        default='$PISA/pisa/images/cfx/',
        help='location onto which to store the plot'
    )
    parser.add_argument(
        '--logv', default=False, action='store_true',
        help='flag to specifiy whether to use a log scale on the z axis'
    )
    parser.add_argument(
        '--ylim', type=float, default=None, nargs=2, metavar=('YMIN', 'YMAX'),
        help='set the limits of the y axis'
    )
    parser.add_argument(
        '-v', '--verbose', action='count', default=0,
        help='set verbosity level'
    )

    args = parser.parse_args()
    return args


def format_2D_axis(ax, binning):
    x_bins, y_bins = binning.dims

    if r'{0:~}'.format(x_bins.units) == '':
        ax.set_xlabel(r'${0}$'.format(x_bins.tex).replace('unsmeared', 'true'),
                      fontsize=18)
    else:
        ax.set_xlabel(r'${0}({1:~})$'.format(
            x_bins.tex, x_bins.units).replace('unsmeared', 'true'), fontsize=18
        )
    if r'{0:~}'.format(y_bins.units) == '':
        ax.set_ylabel(r'${0}$'.format(y_bins.tex).replace('unsmeared', 'true'),
                      fontsize=18)
    else:
        ax.set_ylabel(r'${0}({1:~})$'.format(
            y_bins.tex, y_bins.units).replace('unsmeared', 'true'), fontsize=18
        )

    plt.xticks(x_bins.bin_edges.m, ['%.2f' % x for x in x_bins.bin_edges.m])
    plt.yticks(y_bins.bin_edges.m, ['%.2f' % y for y in y_bins.bin_edges.m])
    ax.tick_params(axis='x', labelsize=10)
    ax.tick_params(axis='y', labelsize=10)

    for xmaj in x_bins.bin_edges.m:
        ax.axvline(x=xmaj, ls='--', color='black', alpha=0.7)
    for ymaj in y_bins.bin_edges.m:
        ax.axhline(y=ymaj, ls='--', color='black', alpha=0.7)

    x_lims = (np.min(x_bins.bin_edges.m), np.max(x_bins.bin_edges.m))
    y_lims = (np.min(y_bins.bin_edges.m), np.max(y_bins.bin_edges.m))
    ax.set_xlim(x_lims)
    ax.set_ylim(y_lims)


def plot_2D(ax, map, colour, logv, vlimits):
    logging.info('Plotting Map {0}'.format(map.name))
    binning = [x.m for x in map.binning.bin_edges]
    hist_T = unp.nominal_values(map.hist).T
    if logv:
        norm = mpl.colors.LogNorm(vmin=vlimits[0], vmax=vlimits[1])
        cax = ax.pcolormesh(
            binning[0], binning[1], hist_T, cmap=plt.get_cmap(colour), norm=norm
        )
    else:
        cax = ax.pcolormesh(
            binning[0], binning[1], hist_T, cmap=plt.get_cmap(colour),
            vmin=vlimits[0], vmax=vlimits[1]
        )

    if map.binning.dims[0].is_log:
        ax.set_xscale('log')
    if map.binning.dims[1].is_log:
        ax.set_yscale('log')

    v = ((vlimits[1] - vlimits[0]) * (6/8.)) + vlimits[0]
    if vlimits[0] == vlimits[1]:
        v = vlimits[0]

    for edex, eedge in enumerate(binning[0][:-1]):
        for zdex, zedge in enumerate(binning[1][:-1]):
            if hist_T[zdex][edex] < v:
                text_color = 'black'
            else:
                text_color = 'white'
            plt.annotate('{0:.3f}'.format(hist_T[zdex][edex]),
                         xy=(binning[0][edex]+(binning[0][edex+1]-binning[0][edex])/8., 
                             binning[1][zdex]+(binning[1][zdex+1]-binning[1][zdex])/3.) ,
                         xycoords='data', color=text_color, fontsize=11)

    return cax


def make_plot(Map, outfile, logv=False, ylim=None):
    fig = plt.figure(figsize=(12, 8))
    fig.suptitle('${0}$'.format(Map.tex), fontsize=14)
    ax = fig.add_subplot(111)

    binning = Map.binning

    vmin = np.min(Map.hist)
    vmax = np.max(Map.hist)
    vlimits = map(unp.nominal_values, (vmin, vmax))

    cax = plot_2D(ax, Map, 'Blues', logv, vlimits)
    format_2D_axis(ax, binning)

    cb = fig.colorbar(cax, ax=ax)
    cb.ax.tick_params(labelsize=11)
    cb.set_label(r'efficiency (%)', fontsize=14)

    fig.savefig(outfile, bbox_inches='tight')

if __name__ == "__main__":
    args = parse_args()
    set_verbosity(args.verbose)

    logging.info('Loading MapSet from file {0}'.format(args.infile))
    input_MapSet = MapSet.from_json(args.infile)

    for Map in input_MapSet:
        outfile = args.outdir + '/' + Map.name
        logging.info('outfile {0}'.format(outfile))
        # turn to percent
        Map *= 100.
        make_plot(
            Map = Map,
            outfile = outfile,
            logv = args.logv,
            ylim = args.ylim,
        )
