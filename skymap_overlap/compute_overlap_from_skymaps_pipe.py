#!/usr/bin/env python
from __future__ import print_function
import os
from pycondor import Job, Dagman
import argparse
import itertools
import shutil

def get_filename_prefix(filename):
    basename = os.path.basename(filename)
    if ".fits.gz" in basename:
        return basename.split(".fits.gz")[0]
    if ".fits" in basename:
        return basename.split(".fits")[0]
    # Fail-safe option
    return os.path.splitext(basename)[0]


def main():
    parser = argparse.ArgumentParser(description = "Compute pair-wise overlap of a batch of skymaps")
    parser.add_argument("--skymap", metavar="PATH", action="append", help="A list of paths pointing to the probability skymaps")
    parser.add_argument("--accounting-tag", type=str, default="ligo.dev.o3.cbc.lensing.multi", help="Accounting tag")
    parser.add_argument("--request-disk", type=float, default=1.0, help="Disk usage (in GB) expected")
    parser.add_argument("--request-memory", type=float, default=4.0, help="Memory usage (in GB) expected")
    parser.add_argument("--slurm", action="store_true", help="Run on a condor+slurm cluster")
    parser.add_argument("--plot", action = "store_true", help = "Visualize the skymaps")
    parser.add_argument("--verbose", action="store_true", help="Be very verbose")

    args = parser.parse_args()

    compute_overlap_job_name = "compute_overlap"
    pairwise_overlap_out_str = "{prefix_1}_{prefix_2}_overlap.dat"

    # Directories for HTCondor
    try:
        os.makedirs(compute_overlap_job_name)
    except:
        pass
    error = os.path.abspath("logs")
    output = os.path.abspath("logs")
    log = os.path.abspath("logs")
    submit = os.path.abspath("")

    # Create a DAG (but actually each node is independent of each other)
    dag = Dagman(
        name="dag_compute_overlap_from_skymaps",
        submit=submit,
    )

    universe = "vanilla"
    extra_lines = [
        "accounting_group = {}".format(args.accounting_tag),
        "request_memory = {:.1f} GB".format(args.request_memory),
        "request_disk = {:.1f} GB".format(args.request_disk),
        "request_cpus = 1",
        "getenv = True",
    ]
    if args.slurm:
        universe = "grid"
        extra_lines.append("grid_resource = batch slurm")

    f = open("compute_overlap_from_skymaps.sh", "w")

    # Compute overlap
    if len(args.skymap) >= 2:
        # At least two skymaps, now we can compute the pairwise overlap
        compute_overlap_job = Job(
            name="job_"+compute_overlap_job_name,
            executable=shutil.which("compute_overlap"),
            universe=universe,
            error=error,
            output=output,
            log=log,
            dag=dag,
            extra_lines=extra_lines,
        )

        for skymap_1, skymap_2 in list(itertools.combinations(args.skymap, 2)):
            prefix_1 = get_filename_prefix(skymap_1)
            prefix_2 = get_filename_prefix(skymap_2)

            argument_str = ""
            if args.verbose:
                argument_str += " --verbose"
            if args.plot:
                argument_str += " --plot"
            argument_str += " --skymap " + os.path.abspath(skymap_1) + " --skymap " + os.path.abspath(skymap_2) + \
                " --output " + os.path.abspath(os.path.join(compute_overlap_job_name, pairwise_overlap_out_str.format(prefix_1=prefix_1, prefix_2=prefix_2)))
            compute_overlap_job.add_arg(argument_str, retry=3)

            f.write(f"compute_overlap {argument_str}\n")

    dag.build(fancyname=False)
    f.close()
