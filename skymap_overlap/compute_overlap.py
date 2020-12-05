#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import healpy as hp
import ligo.skymap.io
import ligo.skymap.postprocess
import ligo.skymap.plot
import matplotlib as mpl
mpl.use("Agg")
from matplotlib import pyplot as plt
from matplotlib import cm as cm
from matplotlib.colors import ListedColormap
import matplotlib.patches as mpatches
from ligo.skymap.bayestar import rasterize
import argparse
import sys
import os

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
    def __init__(self):
        self.name = "cross_hpd_statistic"


    @staticmethod
    def get_ra_dec_from_skymap(skymap):
        index_of_max = np.argmax(skymap)
        nside = hp.npix2nside(len(skymap))
        theta, phi = hp.pix2ang(nside, index_of_max, nest=True)
        return phi,np.pi/2-theta

    def compute_overlap(self,skymap1,skymap2,single_skymap1,single_skymap2):
        from ligo.skymap.postprocess.crossmatch import crossmatch
        from astropy.coordinates import SkyCoord
        ra, dec = self.get_ra_dec_from_skymap(single_skymap1)
        coord = SkyCoord(ra,dec,unit="rad")
        result = crossmatch(skymap2,coord)
        searched_prob_1 = result.searched_prob
        ra, dec = self.get_ra_dec_from_skymap(single_skymap2)
        coord = SkyCoord(ra,dec,unit="rad")
        result = crossmatch(skymap1,coord)
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

def plot_skymaps(skymap_1, skymap_2, label_1="", label_2="", filename="skymaps.pdf"):
    assert len(skymap_1) == len(skymap_2), "The two skymaps should have the same resolution."
    area_per_pixel = hp.nside2pixarea(hp.npix2nside(len(skymap_1)), degrees=True)

    fig = plt.figure()
    ax = plt.axes(projection='astro hours mollweide')
    ax.grid()

    # Plot probability per sq. deg
    skymap_1_persqdeg = skymap_1/area_per_pixel
    skymap_2_persqdeg = skymap_2/area_per_pixel

    cylon = cm.get_cmap('cylon')
    # Make custom colormap
    viridis = cm.get_cmap('viridis')
    viridis_transparent = viridis(np.arange(viridis.N))
    viridis_transparent[:,-1] = np.linspace(0, 1, viridis.N)
    viridis_transparent = ListedColormap(viridis_transparent)

    im_1 = ax.imshow_hpx((skymap_1_persqdeg, 'ICRS'), nested=True, vmin=0., vmax=skymap_1_persqdeg.max(), cmap='cylon')
    im_2 = ax.imshow_hpx((skymap_2_persqdeg, 'ICRS'), nested=True, vmin=0., vmax=skymap_2_persqdeg.max(), cmap=viridis_transparent)

    # Fake the legend
    # Create a patch (proxy artist) for every color
    patches = [mpatches.Patch(color=im_1.cmap(10), label=label_1), mpatches.Patch(color=im_2.cmap(100), label=label_2)]
    plt.legend(handles=patches, bbox_to_anchor=(0.25, 1.2), loc=1, borderaxespad=0.)

    plt.savefig(filename)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Compute overlap in sky given two FITS skymaps")
    parser.add_argument("--skymap", action="append", type=str, metavar="PATH", help="Path to the two sets of FITS skymaps")
    parser.add_argument("--output", type=str, metavar="PATH", help="Path to the text file storing the output")
    parser.add_argument("--plot", action = "store_true", help = "Visualize the skymaps")
    parser.add_argument("--verbose", action = "store_true", help = "Be very verbose")
    args = parser.parse_args()

    if args.skymap is None:
        args.skymap = []
    assert len(args.skymap) == 2, "Currently the script only supports computing the overlap between 2 skymaps. {} given.".format(len(args.skymap))

    posterior_overlap = PosteriorOverlap()
    normalized_posterior_overlap = NormalizedPosteriorOverlap()
    credible_region_overlap = CredibleRegionOverlap(90) # 90% CR overlap
    cross_hpd_statistic = CrossHPDStatistic()

    statistics = [
        posterior_overlap,
        normalized_posterior_overlap,
        credible_region_overlap,
        cross_hpd_statistic,
    ]
    overlap_values = []

    if args.verbose:
        print("Loading {} and {}".format(*args.skymap), file=sys.stdout)

    skymap_1_multi_order = ligo.skymap.io.fits.read_sky_map(args.skymap[0], moc=True)
    skymap_2_multi_order = ligo.skymap.io.fits.read_sky_map(args.skymap[1], moc=True)

    skymap_1, skymap_2 = enforce_same_resolution(
        read_skymap(args.skymap[0]),
        read_skymap(args.skymap[1])
    )

    for stat in statistics:
        if args.verbose:
            print("Calculating {} statistic".format(stat.name), file=sys.stdout)
        if type(stat) == CrossHPDStatistic:
            overlap = stat.compute_overlap(skymap_1_multi_order,skymap_2_multi_order,skymap_1, skymap_2)
        else:
            overlap = stat.compute_overlap(skymap_1, skymap_2)
        overlap_values.append(overlap)

    out_str = ",".join([stat.name for stat in statistics]) + "\n"
    out_str += ",".join(["{}".format(v) for v in overlap_values])

    if args.output is not None:
        with open(args.output, "w") as f:
            f.write(out_str)
    else:
        print(out_str, file=sys.stderr)

    if args.plot:
        skymap_1_label = os.path.basename(args.skymap[0]).split(".fits")[0]
        skymap_2_label = os.path.basename(args.skymap[1]).split(".fits")[0]
        if args.output is not None and args.output.endswith(".dat"):
            out_plot_filename = args.output.replace(".dat", ".pdf")
        else:
            out_plot_filename = "skymaps.pdf"
        plot_skymaps(skymap_1, skymap_2, skymap_1_label, skymap_2_label, out_plot_filename)

main()
