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
    """
    Base class for overlap statistics

    Attributes
    ----------
    name : str

    Methods
    -------
    compute_overlap(*skymaps)
    """
    def __init__(self):
        """
        Constructor
        """
        self.name = "overlap"

    @staticmethod
    def compute_overlap(*skymaps):
        return 0.0


class PosteriorOverlap(OverlapStatistic):
    def __init__(self):
        """
        Constructor
        """
        self.name = "posterior_overlap"

    @staticmethod
    def compute_overlap(*skymaps):
        """
        Compute the overlap between two skymaps

        Parameters
        ----------
        skymaps : list of skymaps

        Returns
        -------
        overlap : overlap
        """
        # Simply sum over all pixels
        # NOTE To avoid under/over-flow, add the log of pdfs then exponentiate
        _out = np.zeros_like(skymaps[0])
        for skymap in skymaps:
            _out += np.log(skymap)
        return np.nansum(np.exp(_out))


class NormalizedPosteriorOverlap(OverlapStatistic):
    """
    Normalized posterior overlap
    """
    def __init__(self):
        """
        Constructor
        """
        self.name = "normalized_posterior_overlap"

    @staticmethod
    def compute_overlap(*skymaps):
        """
        Compute the normalized posterior overlap between two skymaps

        Parameters
        ----------
        skymaps : list of skymaps

        Returns
        -------
        overlap : overlap
        """
        unnormalized_overlap = PosteriorOverlap.compute_overlap(*skymaps)
        normalization = np.prod([np.sqrt(PosteriorOverlap.compute_overlap(m, m)) for m in skymaps])
        return unnormalized_overlap/normalization

class CredibleRegionOverlap(OverlapStatistic):
    def __init__(self, percent):
        """
        Constructor

        Parameters
        ----------
        percent : credible region percentage
        """
        assert 0 < percent < 100, "percent must be between 0 and 100"
        self.name = "{percent}%_credible_region_overlap".format(percent=percent)
        self.percent = percent

    @staticmethod
    def mask_skymap(skymap, percent):
        """
        Mask a skymap with a given credible region percentage

        Parameters
        ----------
        skymap : skymap
        percent : credible region percentage

        Returns
        -------
        masked_skymap : masked skymap
        """
        masked_skymap = np.zeros_like(skymap)
        masked_skymap[ligo.skymap.postprocess.util.find_greedy_credible_levels(skymap) <= percent/100.] = 1.0
        return masked_skymap

    @staticmethod
    def count_masked_pixel(skymap):
        """
        Count the number of masked pixels in a skymap

        Parameters
        ----------
        skymap : skymap

        Returns
        -------
        number of masked pixels
        """
        return len(skymap[skymap == 1.0])

    def compute_overlap(self, *skymaps):
        """
        Compute the overlap between two skymaps

        Parameters
        ----------
        skymaps : list of skymaps

        Returns
        -------
        overlap : overlap
        """
        masked_skymaps = [self.mask_skymap(m, self.percent) for m in skymaps]
        joint_masked_skymaps = np.multiply(*masked_skymaps)
        return self.count_masked_pixel(joint_masked_skymaps)/np.amin([self.count_masked_pixel(m) for m in masked_skymaps])


class CrossHPDStatistic(OverlapStatistic):
    """
    Cross HPD statistic
    """
    def __init__(self):
        """
        Constructor
        """
        self.name = "cross_hpd_statistic"

    @staticmethod
    def get_ra_dec_from_skymap(skymap):
        """
        Get the ra and dec of the maximum probability pixel in a skymap

        Parameters
        ----------
        skymap : skymap

        Returns
        -------
        ra : ra of the maximum probability pixel
        dec : dec of the maximum probability pixel
        """
        index_of_max = np.argmax(skymap)
        nside = hp.npix2nside(len(skymap))
        theta, phi = hp.pix2ang(nside, index_of_max)
        return phi, np.pi/2-theta

    def compute_overlap(self, skymap1, skymap2, single_skymap1, single_skymap2):
        """
        Compute the overlap between two skymaps

        Parameters
        ----------
        skymap1 : skymap
        skymap2 : skymap
        single_skymap1 : single skymap
        single_skymap2 : single skymap

        Returns
        -------
        overlap : overlap
        """
        from ligo.skymap.postprocess.crossmatch import crossmatch
        from astropy.coordinates import SkyCoord
        ra, dec = self.get_ra_dec_from_skymap(single_skymap1)
        coord = SkyCoord(ra, dec, unit="rad")
        result = crossmatch(skymap2, coord)
        searched_prob_1 = np.min([result.searched_prob, 1.0])
        ra, dec = self.get_ra_dec_from_skymap(single_skymap2)
        coord = SkyCoord(ra, dec, unit="rad")
        result = crossmatch(skymap1, coord)
        searched_prob_2 = np.min([result.searched_prob, 1.0])
        return np.max([1-searched_prob_1, 1-searched_prob_2])


def read_skymap(filename):
    """
    Read a skymap from a file

    Parameters
    ----------
    filename : filename

    Returns
    -------
    skymap : skymap
    """
    hpx, _ = ligo.skymap.io.fits.read_sky_map(filename)
    return hpx


def enforce_same_resolution(*skymaps):
    """
    Enforce that all skymaps have the same resolution

    Parameters
    ----------
    skymaps : list of skymaps

    Returns
    -------
    skymaps : list of skymaps
    """
    skymaps = list(skymaps)
    nside_min = np.amin([hp.get_nside(hpx) for hpx in skymaps])
    for i in range(len(skymaps)):
        skymaps[i] = ligo.skymap.postprocess.util.smooth_ud_grade(skymaps[i], nside_min)
    return skymaps


def plot_skymaps(skymaps, labels, cmaps, filename="skymaps.pdf"):
    """
    Plot skymaps

    Parameters
    ----------
    skymaps : list of skymaps
    labels : list of labels
    cmaps : list of cmaps
    filename : filename
    """
    # Sanity check
    assert len(skymaps) == len(labels)
    assert len(skymaps) == len(cmaps)

    # Prepare the custom cmap instances
    bkgrd = np.asarray([1., 1., 1.,])
    cmap_instances = []
    for idx, cmap in enumerate(cmaps):
        cmap = cm.get_cmap(cmap)
        if idx == 0:
            cmap_instances.append(cmap)

        cmap_transparent = cmap(np.arange(cmap.N))
        alphas = np.linspace(0, 1, cmap.N)
        for j in range(cmap.N):
            cmap_transparent[j,:-1] = cmap_transparent[j,:-1]*alphas[j] + bkgrd*(1. - alphas[j])
        cmap_transparent = ListedColormap(cmap_transparent)
        cmap_instances.append(cmap_transparent)

    skymaps = enforce_same_resolution(*skymaps)

    fig = plt.figure()
    ax = plt.axes(projection='astro hours mollweide')
    ax.grid()
    contours = [ax.contour_hpx((ligo.skymap.postprocess.util.find_greedy_credible_levels(skymap), 'ICRS'), nested=False, cmap=cmap_instances[idx]) for idx, skymap in enumerate(skymaps)]
    patches = [mpatches.Patch(color=contours[idx].cmap(10**(idx+1)), label=label) for idx, label in enumerate(labels)]
    for idx, contour in enumerate(contours):
        ax.clabel(contour, inline=True, fontsize=6)
    plt.legend(handles=patches)
    plt.tight_layout()
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
        plot_skymaps([skymap_1, skymap_2], [skymap_1_label, skymap_2_label], ["cylon", "viridis"], out_plot_filename)

