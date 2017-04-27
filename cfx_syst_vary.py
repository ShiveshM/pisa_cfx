#! /usr/bin/env python
# author: S. Mandalia
# date: January 5, 2017
"""
For a pipeline, individually vary each floating systematic by +/- one sigma and
then compare the produced Map to the nominal Map.

"""
import os
import argparse

import numpy as np

from pisa import ureg, Q_
from pisa.core.map import MapSet
from pisa.core.param import Param, ParamSet
from pisa.core.pipeline import Pipeline
from pisa.utils import fileio, resources
from pisa.utils.flavInt import ALL_NUFLAVINTS, NuFlavIntGroup
from pisa.utils.log import logging, set_verbosity

from pisa.scripts.compare import compare

from plot_systematic import make_plot


__all__ = ['SIGMA', 'parse_args', 'main']


SIGMA = -1/2.


def parse_args():
    """Get command line arguments"""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '--outdir',  metavar='DIR', type=str, required=True,
        help='''Store output plots to this directory'''
    )
    parser.add_argument(
        '--cfx-pipeline', type=resources.find_resource, metavar='FILE',
        required=True, help='''Standard CFX pipeline settings config file'''
    )
    parser.add_argument(
        '--signal', type=str, default='numu_cc+numubar_cc',
        help='''Definition of signal with which to unfold'''
    )
    parser.add_argument(
        '--sigma', type=float, default=1., metavar='FLOAT',
        help='''Size of systematic shift in units sigma'''
    )
    parser.add_argument(
        '--json', action='store_true',
        help='''Save output maps in compressed json (json.bz2) format.'''
    )
    parser.add_argument(
        '--pdf', action='store_true',
        help='''Save plots in PDF format. If neither this nor --png is
        specified, no plots are produced.'''
    )
    parser.add_argument(
        '--png', action='store_true',
        help='''Save plots in PNG format. If neither this nor --pdf is
        specfied, no plots are produced.'''
    )
    parser.add_argument(
        '--center-zero', action='store_true',
        help='center the z axis on 0 when plotting'
    )
    parser.add_argument(
        '--diff-min', type=float, required=False,
        help='''Difference plot vmin; if you specify only one of --diff-min or
        --diff-max, symmetric limits are automatically used (min = -max).'''
    )
    parser.add_argument(
        '--diff-max', type=float, required=False,
        help='''Difference plot max; if you specify only one of --diff-min or
        --diff-max, symmetric limits are automatically used (min = -max).'''
    )
    parser.add_argument(
        '--fract-diff-min', type=float, required=False,
        help='''Fractional difference plot vmin; if you specify only one of
        --fract-diff-min or --fract-diff-max, symmetric limits are
        automatically used (min = -max).'''
    )
    parser.add_argument(
        '--fract-diff-max', type=float, required=False,
        help='''Fractional difference plot max; if you specify only one of
        --fract-diff-min or --fract-diff-max, symmetric limits are
        automatically used (min = -max).'''
    )
    parser.add_argument(
        '--asymm-min', type=float, required=False,
        help='''Asymmetry plot vmin; if you specify only one of --asymm-min or
        --asymm-max, symmetric limits are automatically used (min = -max).'''
    )
    parser.add_argument(
        '--asymm-max', type=float, required=False,
        help='''Fractional difference plot max; if you specify only one of
        --asymm-min or --asymm-max, symmetric limits are automatically used
        (min = -max).'''
    )
    parser.add_argument(
        '-v', action='count',
        help='Set verbosity level; repeat -v for higher level.'
    )
    args = parser.parse_args()
    return args


def main():
    global SIGMA
    args = vars(parse_args())
    set_verbosity(args.pop('v'))
    center_zero = args.pop('center_zero')

    make_pdf = False
    if args['pdf']:
        make_pdf = True
        args['pdf']= False

    outdir = args.pop('outdir')
    fileio.mkdir(outdir, mode=0755)
    SIGMA *= args.pop('sigma')

    cfx_pipe = Pipeline(args.pop('cfx_pipeline'))

    signal = args.pop('signal').replace(' ', '').split(',')
    output_str = []
    for name in signal:
        if 'muons' in name or 'noise' in name:
            raise AssertionError('Are you trying to unfold muons/noise?')
        elif 'all_nu' in name:
            output_str = [str(NuFlavIntGroup(f)) for f in ALL_NUFLAVINTS]
        else:
            output_str.append(NuFlavIntGroup(name))
    output_str = [str(f) for f in output_str]
    cfx_pipe._output_names = output_str

    # Turn off stat fluctuations
    stat_param = cfx_pipe.params['stat_fluctuations']
    stat_param.value = 0 * ureg.dimensionless
    cfx_pipe.update_params(stat_param)

    # Get nominal Map
    re_param = cfx_pipe.params['regularisation']
    re_param.value = 0 * ureg.dimensionless
    cfx_pipe.update_params(re_param)
    nom_out = cfx_pipe.get_outputs()

    re_param.reset()
    cfx_pipe.update_params(re_param)

    params = ParamSet()
    for param in cfx_pipe.params:
        if param.name != 'dataset':
            params.extend(param)

    free = params.free
    logging.info('Free params = {0}'.format(free))
    contin = True
    for f in free:
        if 'hole_ice' not in f.name and 'dom_eff' not in f.name:
            continue
        # if 'atm_muon_scale' in f.name:
        #     contin = False
        # if contin:
        #     continue

        logging.info('Working on parameter {0}'.format(f.name))
        if f.prior.kind != 'uniform':
            # Use deltaLLH = SIGMA to define +/- sigma for non-uniform
            scan_over = np.linspace(*f.range, num=1000) * f.range[0].u
            llh = f.prior.llh(scan_over)
            dllh = llh - np.min(-llh)

            mllh_idx = np.argmin(-llh)
            if mllh_idx == 0:
                l_sig_idx = 0
            else:
                l_sig_idx = np.argmin(np.abs(dllh[:mllh_idx] - SIGMA))
            u_sig_idx = np.argmin(np.abs(dllh[mllh_idx:] - SIGMA)) + mllh_idx

            l_sigma = scan_over[l_sig_idx]
            u_sigma = scan_over[u_sig_idx]
        else:
            l_sigma = f.range[0]
            u_sigma = f.range[1]

        logging.info('Setting {0} lower sigma bound to '
                     '{1}'.format(f.name, l_sigma))
        f.value = l_sigma
        cfx_pipe.update_params(f)
        l_out = cfx_pipe.get_outputs()

        logging.info('Setting {0} upper sigma bound to '
                     '{1}'.format(f.name, u_sigma))
        f.value = u_sigma
        cfx_pipe.update_params(f)
        u_out = cfx_pipe.get_outputs()

        f.reset()
        cfx_pipe.update_params(f)

        f_outdir = outdir+'/'+f.name
        l_outdir = f_outdir+'/'+'lower'
        u_outdir = f_outdir+'/'+'upper'
        fileio.mkdir(f_outdir)
        fileio.mkdir(l_outdir)
        fileio.mkdir(u_outdir)

        compare(
            outdir=l_outdir,
            ref=MapSet([nom_out]),
            ref_label='baseline',
            test=MapSet([l_out]),
            test_label=r'-sigma',
            **args
        )
        compare(
            outdir=u_outdir,
            ref=MapSet([nom_out]),
            ref_label='baseline',
            test=MapSet([u_out]),
            test_label=r'+sigma',
            **args
        )

        l_in_mapset = l_outdir+'/'+'fract_diff__-sigma___baseline.json.bz2'
        u_in_mapset = u_outdir+'/'+'fract_diff__+sigma___baseline.json.bz2'
        l_in_map = MapSet.from_json(l_in_mapset).pop() * 100.
        u_in_map = MapSet.from_json(u_in_mapset).pop() * 100.

        if make_pdf:
            outfile = f_outdir + '/systematic_effect.pdf'
        else:
            outfile = f_outdir + '/systematic_effect.png'
        title = r'% effect on ' + r'${0}$'.format(l_in_map.tex) + \
                ' event counts for {0} parameter'.format(f.name)
        sub_titles = (
            r'(-\sigma - {\rm baseline}) \:/\: {\rm baseline}',
            r'(+\sigma - {\rm baseline}) \:/\: {\rm baseline}'
        )
        make_plot(
            maps=(l_in_map, u_in_map),
            outfile=outfile,
            logv=False,
            center_zero=center_zero,
            vlabel=r'({\rm change} - {\rm baseline}) \:/\: {\rm baseline} (%)',
            title=title,
            sub_titles=sub_titles
        )


main.__doc__ = __doc__


if __name__ == '__main__':
    main()
