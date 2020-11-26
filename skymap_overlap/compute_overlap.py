#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import healpy as hp
import ligo.skymap.io
import ligo.skymap.postprocess
from ligo.skymap.postprocess.crossmatch import crossmatch
from astropy.coordinates import SkyCoord
import argparse
import sys
import os
import bilby

class OverlapStatistic(object):
    def __init__(self):
        self.name = "overlap"

    @staticmethod
    def compute_overlap(*skymaps):
        return 0.0

class PosteriorOverlap(OverlapStatistic):
    def __init__(self):
        self.name = "posterior_overlap"

    @staticmethod
    def compute_overlap(*skymaps):
        # Simply sum over all pixels
        return np.sum(np.multiply(*skymaps))

class NormalizedPosteriorOverlap(OverlapStatistic):
    def __init__(self):
        self.name = "normalized_posterior_overlap"

    @staticmethod
    def compute_overlap(*skymaps):
        unnormalized_overlap = PosteriorOverlap.compute_overlap(*skymaps)
        normalization = np.prod([np.sqrt(PosteriorOverlap.compute_overlap(m, m)) for m in skymaps])
        return unnormalized_overlap/normalization

class CredibleRegionOverlap(OverlapStatistic):
    def __init__(self, percent):
        assert 0 < percent < 100, "percent must be between 0 and 100"
        self.name = "{percent}%_credible_region_overlap".format(percent=percent)
        self.percent = percent

    @staticmethod
    def mask_skymap(skymap, percent):
        masked_skymap = np.zeros_like(skymap)
        masked_skymap[ligo.skymap.postprocess.util.find_greedy_credible_levels(skymap) <= percent/100.] = 1.0
        return masked_skymap

    @staticmethod
    def count_masked_pixel(skymap):
        return len(skymap[skymap == 1.0])

    def compute_overlap(self, *skymaps):
        masked_skymaps = [self.mask_skymap(m, self.percent) for m in skymaps]
        joint_masked_skymaps = np.multiply(*masked_skymaps)
        return self.count_masked_pixel(joint_masked_skymaps)/np.amin([self.count_masked_pixel(m) for m in masked_skymaps])

class CrossHPDStatistic(OverlapStatistic):
    def __init__(self, skymap_path1, skymap_path2):
        self.name = "cross_hpd_statistic"
        self.path1 = os.path.dirname(skymap_path1)
        self.json_name1 = os.path.basename(skymap_path1).split(".fits")[0]
        self.path2 = os.path.dirname(skymap_path2)
        self.json_name2 = os.path.basename(skymap_path2).split(".fits")[0]
        self.skymap1 = ligo.skymap.io.fits.read_sky_map(skymap_path1,moc=True)
        self.skymap2 = ligo.skymap.io.fits.read_sky_map(skymap_path2,moc=True)

    @staticmethod
    def get_ra_dec_from_json(path, name):
        result = bilby.result.read_in_result(path+'/'+name+'.json')
        result.posterior['log_posterior'] = result.posterior['log_likelihood'] + result.posterior['log_prior']
        ra,dec = result.posterior.sort_values(by="log_posterior").iloc[-1][['ra','dec']]
        return ra,dec

    def compute_overlap(self,*skymaps):
        ra, dec = self.get_ra_dec_from_json(self.path1,self.json_name1)
        coord = SkyCoord(ra,dec,unit="rad")
        result = crossmatch(self.skymap1,coord)
        searched_prob_1 = result.searched_prob
        ra, dec = self.get_ra_dec_from_json(self.path2,self.json_name2)
        coord = SkyCoord(ra,dec,unit="rad")
        result = crossmatch(self.skymap2,coord)
        searched_prob_2 = result.searched_prob
        return np.max([1-searched_prob_1, 1-searched_prob_2])




def read_skymap(filename):
    hpx, _ = ligo.skymap.io.fits.read_sky_map(filename)
    return hpx

def enforce_same_resolution(*skymaps):
    skymaps = list(skymaps)
    nside_min = np.amin([hp.get_nside(hpx) for hpx in skymaps])
    for i in range(len(skymaps)):
        skymaps[i] = ligo.skymap.postprocess.util.smooth_ud_grade(skymaps[i], nside_min)
        
    return skymaps

def main():
    parser = argparse.ArgumentParser(description="Compute overlap in sky given two FITS skymaps")
    parser.add_argument("--skymap", action="append", type=str, metavar="PATH", help="Path to the two sets of FITS skymaps")
    parser.add_argument("--output", type=str, metavar="PATH", help="Path to the text file storing the output")
    parser.add_argument("--verbose", action = "store_true", help = "Be very verbose")
    args = parser.parse_args()

    if args.skymap is None:
        args.skymap = []
    assert len(args.skymap) == 2, "Currently the script only supports computing the overlap between 2 skymaps. {} given.".format(len(args.skymap))

    posterior_overlap = PosteriorOverlap()
    normalized_posterior_overlap = NormalizedPosteriorOverlap()
    credible_region_overlap = CredibleRegionOverlap(90) # 90% CR overlap
    cross_hpd_statistic = CrossHPDStatistic(args.skymap[0],args.skymap[1])

    statistics = [
        posterior_overlap,
        normalized_posterior_overlap,
        credible_region_overlap,
        cross_hpd_statistic,
    ]
    overlap_values = []

    if args.verbose:
        print("Loading {} and {}".format(*args.skymap), file=sys.stdout)

    skymap_1, skymap_2 = enforce_same_resolution(
        read_skymap(args.skymap[0]),
        read_skymap(args.skymap[1])
    )

    for stat in statistics:
        if args.verbose:
            print("Calculating {} statistic".format(stat.name), file=sys.stdout)
        overlap = stat.compute_overlap(skymap_1, skymap_2)
        overlap_values.append(overlap)

    out_str = ",".join([stat.name for stat in statistics]) + "\n"
    out_str += ",".join(["{}".format(v) for v in overlap_values])

    if args.output is not None:
        with open(args.output, "w") as f:
            f.write(out_str)
    else:
        print(out_str, file=sys.stderr)

main()
