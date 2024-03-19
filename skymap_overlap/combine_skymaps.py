#!/usr/bin/env python

# NOTE This is modified by ligo-skymap-combine

from ligo.skymap.tool import ArgumentParser, FileType


def parser():
    parser = ArgumentParser()
    parser.add_argument('input', metavar='INPUT.fits[.gz]',
                        type=FileType('rb'), nargs='+',
                        help='Input sky localizations')
    # FIXME the output option has type str because astropy.io.fits.writeto()
    # only honors the .gz extension when given a file name string (as of 3.0.1)
    parser.add_argument('output', metavar='OUTPUT.fits[.gz]', type=str,
                        help='Output combined sky localization')
    parser.add_argument('--origin', type=str,
                        help='Optional tag describing the organization'
                             ' responsible for the combined output')
    return parser


def main(args=None):
    args = parser().parse_args(args)

    from textwrap import wrap
    import numpy as np
    import astropy_healpix as ah
    from astropy.io import fits
    from astropy.time import Time
    import healpy as hp

    from ligo.skymap.io import read_sky_map, write_sky_map

    input_skymaps = []
    for input_file in args.input:
        with fits.open(input_file) as hdus:
            header = hdus[0].header.copy()
            header.extend(hdus[1].header)

            # Explicitly disable distance map
            data, meta = read_sky_map(hdus, nest=True,
                                      distances=False)

            data = (data,)

        nside = ah.npix_to_nside(len(data[0]))
        input_skymaps.append((nside, data[0], meta, header))

    max_nside = max(x[0] for x in input_skymaps)

    # upsample sky posteriors to maximum resolution and combine them
    combined_prob = None
    for nside, prob, _, _ in input_skymaps:
        if nside < max_nside:
            prob = hp.ud_grade(prob, max_nside, order_in='NESTED',
                               order_out='NESTED')
        if combined_prob is None:
            combined_prob = np.ones_like(prob)
        combined_prob *= prob

    # normalize joint posterior
    norm = combined_prob.sum()
    if norm == 0:
        raise RuntimeError('input sky localizations are disjoint')
    combined_prob /= norm

    out_kwargs = {'gps_creation_time': Time.now().gps,
                  'nest': True}
    if args.origin is not None:
        out_kwargs['origin'] = args.origin

    out_data = combined_prob

    # save input headers in output history
    out_kwargs['HISTORY'] = []
    for i, x in enumerate(input_skymaps):
        out_kwargs['HISTORY'].append('')
        out_kwargs['HISTORY'].append(
            'Headers of HDUs 0 and 1 of input file {:d}:'.format(i))
        out_kwargs['HISTORY'].append('')
        for line in x[3].tostring(sep='\n',
                                  endcard=False,
                                  padding=False).split('\n'):
            out_kwargs['HISTORY'].extend(wrap(line, 72))

    write_sky_map(args.output, out_data, **out_kwargs)