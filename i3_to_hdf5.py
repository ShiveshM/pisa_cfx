#!/usr/bin/env python
"""
Quick script to convert i3 file to hdf5
"""

import os
import argparse
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

import numpy as np
import h5py

import icecube
from icecube import dataclasses, dataio


class FullPaths(argparse.Action):
    """
    Append user- and relative-paths
    """
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))


def is_dir(dirname):
    """
    Checks if a path is an actual directory, if not found then it creates the
    directory
    """
    if not os.path.isdir(dirname):
        fileio.mkdir(dirname)
    return dirname


def parse_args():
    """Get command line arguments"""
    parser = ArgumentParser(
        description='''Takes an outputted json file containing a MapSet
        object with a single Map in it and makes a plot of it.''',
        formatter_class=ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        '-i', '--indir', type=is_dir, action=FullPaths, metavar='DIR',
        required=True, help='location containing the I3 files'
    )

    parser.add_argument(
        '-o', '--outname', metavar='FILE', type=str, required=True,
        help='output hdf5 filename'
    )

    parser.add_argument(
        '--nu_frac', type=float, default=0.7,
        help='fraction of events that are nu - rest are taken to be nubar'
    )
    args = parser.parse_args()
    return args


def load_file(path, s_parms, nu_frac=1.):
    """
    Get nu parameters from an I3 file
    """
    print 'Loading I3 file {0}'.format(path)
    infile = dataio.I3File(path, 'r')

    parms = {}
    for p in s_parms:
        parms[p] = []

    while(infile.more()):
        frame = infile.pop_daq()
        if str(frame) == 'None':
            infile.rewind()
            break

        mctree = frame['I3MCTree']
        nu = dataclasses.get_most_energetic_neutrino(mctree)
        energy = nu.energy
        coszen = np.cos(nu.dir.zenith)
        azimuth = nu.dir.azimuth * 180. / np.pi
        ptype = nu.type.real

        ow = frame['I3MCWeightDict']['OneWeight']
        nevents = frame['I3MCWeightDict']['NEvents']
        oneweight = ow / nevents
        if ptype > 0:
            oneweight /= nu_frac
        else:
            oneweight /= (1 - nu_frac)

        interaction = frame['I3MCWeightDict']['InteractionType']
        volume = frame['I3MCWeightDict']['GeneratorVolume']

        parms['energy'].append(energy)
        parms['coszen'].append(coszen)
        parms['azimuth'].append(azimuth)
        parms['ptype'].append(ptype)
        parms['oneweight'].append(oneweight)
        parms['interaction'].append(interaction)
        parms['volume'].append(volume)
    return parms

def main():
    args = parse_args()

    s_parms = ['energy', 'coszen', 'azimuth', 'ptype', 'oneweight',
             'interaction', 'volume']
    parms = {}
    for p in s_parms:
        parms[p] = []

    print 'Loading I3 files from dir {0}'.format(args.indir)
    n_files = 0
    for file_name in sorted(os.listdir(args.indir)):
        file_path = args.indir + '/' + file_name
        if os.path.isfile(file_path) \
           and (file_path.endswith('.i3.bz2') or file_path.endswith('.i3.gz')):
            n_files += 1
            sub_parms = load_file(file_path, s_parms, nu_frac=args.nu_frac)
            for key in parms.iterkeys():
                parms[key] += sub_parms[key]

    for p in parms.iterkeys():
        parms[p] = np.array(parms[p])
    parms['oneweight'] /= float(n_files)
    print 'Done loading params from I3 file'

    print 'Saving in file {0}'.format(args.outname)
    h5file = h5py.File(args.outname, 'w')
    for p in parms.iterkeys():
        h5file.create_dataset(p, data=parms[p])
    h5file.close()


if __name__ == "__main__":
    main()
