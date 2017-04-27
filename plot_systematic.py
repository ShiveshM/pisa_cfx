#!/usr/bin/env python
"""
From 2 MapSet json files, create two 2D plots, side by side.
"""

import os
import argparse
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

import numpy as np
from uncertainties import unumpy as unp
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
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
        '-i', '--infiles', type=str, action='append',
        metavar=('FILE', 'FILE'), required=True, nargs=2,
        help='''Location/filenames of the json files containing the MapSet.'''
    )
    parser.add_argument(
        '-o', '--outdir', type=fileio.is_dir, action=FullPaths, metavar='DIR',
        required=True, help='location onto which to store the plot'
    )
    parser.add_argument(
        '-n', '--outname', metavar='FILE', type=str,
        default='untitled.png', help='output filename'
    )
    parser.add_argument(
        '--logv', default=False, action='store_true',
        help='flag to specifiy whether to use a log scale on the y axis'
    )
    parser.add_argument(
        '--center-zero', action='store_true',
        help='center the z axis on 0'
    )
    parser.add_argument(
        '--vlabel', type=str, default=None, metavar='STR',
        help='set the label of the z axis'
    )
    parser.add_argument(
        '--title', type=str, default=None, metavar='STR', help='Title of plot'
    )
    parser.add_argument(
        '--subtitles', type=str, default=None, nargs=2, metavar=('STR', 'STR'),
        help='Title of each subplot'
    )
    parser.add_argument(
        '-v', '--verbose', action='count', default=0,
        help='set verbosity level'
    )
    args = parser.parse_args()
    return args


def format_2D_axis(ax, binning, title):
    x_bins, y_bins = binning.dims

    if r'{0:~}'.format(x_bins.units) == '':
        ax.set_xlabel(r'${0}$'.format(x_bins.tex), fontsize=18)
    else:
        ax.set_xlabel(r'${0}({1:~})$'.format(
            x_bins.tex, x_bins.units), fontsize=18
        )
    if r'{0:~}'.format(y_bins.units) == '':
        ax.set_ylabel(r'${0}$'.format(y_bins.tex), fontsize=18)
    else:
        ax.set_ylabel(r'${0}({1:~})$'.format(
            y_bins.tex, y_bins.units), fontsize=18
        )
    ax.set_title(r'${0}$'.format(title), fontsize=18)

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

    v = ((vlimits[1] - vlimits[0]) * (7/8.)) + vlimits[0]
    if vlimits[0] == vlimits[1]:
        v = vlimits[0]
    f = ((vlimits[1] - vlimits[0]) * (1/8.)) + vlimits[0]
    if vlimits[0] == vlimits[1]:
        f = vlimits[0]

    for edex, eedge in enumerate(binning[0][:-1]):
        for zdex, zedge in enumerate(binning[1][:-1]):
            if hist_T[zdex][edex] > f and hist_T[zdex][edex] < v:
                text_color = 'black'
            else:
                text_color = 'white'
            plt.annotate('{0:.2f}'.format(hist_T[zdex][edex]),
                         xy=(binning[0][edex]+(binning[0][edex+1]-binning[0][edex])/8., 
                             binning[1][zdex]+(binning[1][zdex+1]-binning[1][zdex])/3.) ,
                         xycoords='data', color=text_color, fontsize=11)

    return cax


def make_plot(maps, outfile, logv=False, center_zero=False, vlabel=None,
              title=None, sub_titles=None):
    fig = plt.figure(figsize=(24, 8))
    if title is not None:
        fig.suptitle(title, fontsize=14)
    gs = gridspec.GridSpec(1, 2, width_ratios=[8, 10])
    gs.update(hspace=0.4, wspace=0.15)

    binning = maps[0].binning
    assert binning == maps[1].binning

    vmin = np.min([np.min(maps[0].hist), np.min(maps[1].hist)])
    vmax = np.max([np.max(maps[0].hist), np.max(maps[1].hist)])
    if center_zero:
        max_diff = unp.nominal_values(np.max([vmin, vmax]))
        vlimits = [-max_diff, max_diff]
    else:
        vlimits = map(unp.nominal_values, (vmin, vmax))

    ax0 = fig.add_subplot(gs[0])
    cax = plot_2D(ax0, maps[0], 'RdYlBu', logv, vlimits)
    format_2D_axis(ax0, binning, sub_titles[0])

    ax1 = fig.add_subplot(gs[1])
    cax = plot_2D(ax1, maps[1], 'RdYlBu', logv, vlimits)
    format_2D_axis(ax1, binning, sub_titles[1])

    cb = fig.colorbar(cax, ax=ax1)
    cb.ax.tick_params(labelsize=11)
    if vlabel is not None:
        cb.set_label(r'${0}$'.format(vlabel), fontsize=14)

    fig.savefig(outfile, bbox_inches='tight')


def main():
    args = parse_args()
    set_verbosity(args.verbose)

    assert len(args.infiles[0]) == 2
    assert len(args.subtitles) == 2
    logging.info('Loading MapSet from files {0}'.format(args.infiles[0]))
    o_mapset = MapSet.from_json(args.infiles[0][0])
    t_mapset = MapSet.from_json(args.infiles[0][1])
    assert len(o_mapset) == 1 and len(t_mapset) == 1
    o_map = o_mapset.pop()
    t_map = t_mapset.pop()

    outfile = args.outdir + '/' + args.outname
    make_plot(
        maps = (o_map, t_map),
        outfile = outfile,
        logv = args.logv,
        center_zero = args.center_zero,
        vlabel = args.vlabel,
        title = args.title,
        sub_titles = args.subtitles,
    )


main.__doc__ = __doc__


if __name__ == '__main__':
    main()
