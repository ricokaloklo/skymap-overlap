#!/usr/bin/env python
import numpy as np
import argparse

def compute_overlap(ra_samples_1, dec_samples_1, ra_samples_2, dec_samples_2, nbins=500):
    """
    Compute the overlap between skymap 1 and skymap 2 given the (RA, DEC) samples generated from
    the two probability skymaps p(\Omega) and q(\Omega)
    by \int_{whole sky} p(RA, DEC) q(RA, DEC) d \Omega = \int_{whole sky} p(RA, DEC) q(RA, DEC) d(RA) d(sin DEC)
    """
    assert nbins > 1, "Number of bins must be greater than 1"
    ra_grid = np.linspace(0, 2.0*np.pi, num=nbins)
    sin_dec_grid = np.linspace(-1, 1, num=nbins)
    d_ra = (ra_grid[-1] - ra_grid[0])/nbins
    d_sin_dec = (sin_dec_grid[-1] - sin_dec_grid[0])/nbins

    pdf_1, xedges_1, yedges_1 = np.histogram2d(ra_samples_1, np.sin(dec_samples_1), bins=(ra_grid, sin_dec_grid), normed=True)
    pdf_2, xedges_2, yedges_2 = np.histogram2d(ra_samples_2, np.sin(dec_samples_2), bins=(ra_grid, sin_dec_grid), normed=True)

    return np.sum(np.sum(pdf_1*pdf_2))*d_ra*d_sin_dec


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute overlap in sky given posterior samples of (RA, DEC)")
    parser.add_argument("posterior_samples", nargs=2, type=str, metavar="PATH", help="Path to the two sets of posterior samples of (RA, DEC)")
    parser.add_argument("--output", type=str, metavar="PATH", help="Path to the text file storing the output")
    args = parser.parse_args()

    # Load posterior samples from files
    dec_1, ra_1 = np.transpose(np.loadtxt(args.posterior_samples[0]))
    dec_2, ra_2 = np.transpose(np.loadtxt(args.posterior_samples[1]))

    # Actually compute the overlap
    overlap = compute_overlap(ra_1, dec_1, ra_2, dec_2)

    # Save the output (just a number really)
    np.savetxt(args.output, overlap)