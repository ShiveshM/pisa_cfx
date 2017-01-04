from copy import deepcopy

import numpy as np
import pint
from uncertainties import ufloat
from uncertainties import unumpy as unp

from pisa import ureg, Q_
from pisa.core.param import Param
from pisa.core.pipeline import Pipeline
from pisa.core.map import Map, MapSet
from pisa.utils.log import set_verbosity


if __name__ == '__main__':
    outname = 'test'

    if 'test' in outname:
        set_verbosity(0)
    else:
        set_verbosity(1)

    # livetimes = [2, 3, 4, 5, 6, 7, 8] * ureg.common_year
    livetimes = [1, 4, 16, 64] * ureg.common_year
    # livetimes = [1, 2] * ureg.common_year

    # TODO(shivesh): find a way to sync these automatically
    pipeline = Pipeline('settings/pipeline/example_cfx.cfg')
    gen_pipe = Pipeline('settings/pipeline/gen_lvl.cfg')

    # # sync dataset
    # p_dsc = pipeline.params['data_sample_config']
    # g_dsc = gen_pipe.params['data_sample_config']
    # g_dsc = p_dsc
    # gen_pipe.update_params(g_dsc)

    re_param = pipeline.params['regularisation']
    sf_param = pipeline.params['stat_fluctuations']
    lt_param = pipeline.params['livetime']

    unfold_pipeline_cfg = pipeline.params['unfold_pipeline_cfg']

    mean_perpe = []
    mean_perbin = []
    for idx, lt in enumerate(livetimes):
        print '==========='
        print 'livetime = {0}'.format(lt)
        print '==========='
        mean_perpe.append([])

        lt_param.value = lt
        pipeline.update_params(lt_param)
        gen_pipe.update_params(lt_param)

        u_pipe = Param(name='unfold_pipeline_cfg', value=gen_pipe,
                       is_fixed=True, prior=None, range=None)
        unfold_pipeline_cfg = u_pipe
        pipeline.update_params(unfold_pipeline_cfg)

        # Get nominal
        re_param.value = 0 * ureg.dimensionless
        pipeline.update_params(re_param)
        nom_out = pipeline.get_outputs().pop()

        re_param.value = 2 * ureg.dimensionless
        sf_param.value = 1234 * ureg.dimensionless
        pipeline.update_params(re_param)
        pipeline.update_params(sf_param)

        fe = []
        if 'test' in outname:
            n_trials = 5
        else:
            n_trials = 200
        for x in xrange(n_trials):
            temp_out = pipeline.get_outputs().pop()
            nan_mask = (nom_out.hist == 0)
            div = temp_out.hist[~nan_mask] / nom_out.hist[~nan_mask]
            fe.append(div)
        for f in fe:
            mean_perpe[idx].append(np.mean(f))
            # t = ufloat(np.mean(unp.nominal_values(f)),
            #            np.std(unp.nominal_values(f)))
            mean_perpe[idx].append(t)
        mean_perbin.append(np.mean(fe, axis=0).flatten())
        # a = unp.uarray(np.mean(unp.nominal_values(fe), axis=0).flatten(),
        #                np.std(unp.nominal_values(fe), axis=0).flatten())
        mean_perbin.append(a)

    fe = zip(*mean_perpe)
    print fe

    import matplotlib as mpl
    # headless mode
    mpl.use('Agg')
    # fonts
    mpl.rcParams['mathtext.fontset'] = 'custom'
    mpl.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
    mpl.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
    mpl.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
    from matplotlib import pyplot as plt
    from matplotlib.offsetbox import AnchoredText

    binning = livetimes.m
    fig = plt.figure(figsize=(9, 5))
    ax = fig.add_subplot(111)
    ax.set_xlim(np.min(binning)-1, np.max(binning)+1)
    # ax.set_ylim(0.9, 1.1)
    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=12)

    ax.set_xlabel('livetime (years)', fontsize=18)
    ax.set_ylabel('mean ratio per pe unfolded vs. truth (200 trials)', fontsize=15)

    def get_edges_from_cen(bincen):
        hwidth = 0.5*(bincen[1] - bincen[0])
        return np.append([bincen[0]-hwidth], bincen[:]+hwidth)

    for f in fe:
        fe_0 = np.concatenate([[f[0]], f])
        ax.errorbar(
            binning, unp.nominal_values(f), xerr=0,
            yerr=unp.std_devs(f), capsize=3, alpha=0.5, linestyle='--',
            markersize=2, linewidth=1
        )

    for ymaj in ax.yaxis.get_majorticklocs():
        ax.axhline(y=ymaj, ls=':', color='gray', alpha=0.7, linewidth=1)
    for xmaj in ax.xaxis.get_majorticklocs():
        ax.axvline(x=xmaj, ls=':', color='gray', alpha=0.7, linewidth=1)
    fig.savefig('./images/perpe/'+outname+'.png', bbox_inches='tight', dpi=150)

    mean = zip(*mean_perbin)
    print np.array(mean)

    fig_2 = plt.figure(figsize=(9, 5))
    ax_2 = fig_2.add_subplot(111)
    ax_2.set_xlim(np.min(binning)-1, np.max(binning)+1)
    # ax_2.set_ylim(0.5, 1.5)
    ax_2.tick_params(axis='x', labelsize=14)
    ax_2.tick_params(axis='y', labelsize=12)

    ax_2.set_xlabel('livetime (years)', fontsize=18)
    ax_2.set_ylabel('mean ratio per bin unfolded vs. truth (200 trials)', fontsize=15)

    for f in mean:
        ax_2.errorbar(
            binning, unp.nominal_values(f), xerr=0,
            yerr=unp.std_devs(f), capsize=3, alpha=0.5, linestyle='--',
            markersize=2, linewidth=1
        )

    for ymaj in ax_2.yaxis.get_majorticklocs():
        ax_2.axhline(y=ymaj, ls=':', color='gray', alpha=0.7, linewidth=1)
    for xmaj in ax_2.xaxis.get_majorticklocs():
        ax_2.axvline(x=xmaj, ls=':', color='gray', alpha=0.7, linewidth=1)
    fig_2.savefig('./images/perbin/'+outname+'.png', bbox_inches='tight', dpi=150)
