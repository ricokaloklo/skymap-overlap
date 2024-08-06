#!/usr/bin/env python
from __future__ import print_function
from ligo.gracedb.rest import GraceDb
import argparse
import sys

def download_skymap(id, db, args, use_bayestar_only=False):
	filename = None
	try:
		# Use PE skymap if possible, except when use_bayestar_only is set
		if not use_bayestar_only:
			try:
				# Try bilby first
				filename = "Bilby.multiorder.fits"
				r = db.files(id, filename)
				if args.verbose:
					print("Using Bilby skymap for {0}".format(id), file=sys.stderr)
			except:
				# Fallback to LALInference for older events
				filename = "LALInference.fits.gz"
				r = db.files(id, filename)
				if args.verbose:
					print("Using LALInference skymap for {0}".format(id), file=sys.stderr)
		else:
			filename = ""
			raise TypeError
	except:
		try:
			filename = "bayestar.fits.gz"
			r = db.files(id, filename)
			if args.verbose:
				print("Using Bayestar skymap for {0}".format(id), file=sys.stderr)
		except:
			filename = "subthreshold.bayestar.fits.gz"
			r = db.files(id, filename)
			if args.verbose:
				print("Using subthreshold Bayestar skymap for {0}".format(id), file=sys.stderr)

	out_filename = "{0}_skymap.fits.gz".format(id)
	if filename.endswith("fits.gz"):
		outfile = open(out_filename, "wb")
	else:
		import gzip
		outfile = gzip.open(out_filename, "wb")
	if args.verbose:
		print("Saving skymap to {0}".format(out_filename), file=sys.stderr)
	outfile.write(r.read())
	outfile.close()

def main():
	parser = argparse.ArgumentParser(description = "Download skymaps from a list of events")
	parser.add_argument("event", nargs="+", help = "A list of gravitational-wave events, can be either GID for GW event or SID for superevent")
	parser.add_argument("--bayestar", action="store_true", help="Use bayestar skymap only")
	parser.add_argument("--verbose", action = "store_true", help = "Be very verbose")
	args = parser.parse_args()
	# FIXME Make sure that you have a valid proxy
	client = GraceDb()

	for event_id in args.event:
		try:
			download_skymap(event_id, client, args, use_bayestar_only=args.bayestar)
		except:
			if args.verbose:
				print("Failed to download the skymap for {}".format(event_id), file=sys.stderr)
