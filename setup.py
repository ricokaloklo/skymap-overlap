import setuptools

# Read hanabi/_version.py
# Code from StackOverflow: https://stackoverflow.com/questions/458550/standard-way-to-embed-version-into-python-package
import re
VERSIONFILE="skymap_overlap/_version.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))

setuptools.setup(
    name="skymap-overlap",
    version=verstr,
    author="Rico Ka Lok Lo",
    author_email="kllo@caltech.edu",
    description="Compute overlap between two skymaps",
    long_description="LONG DESCRIPTION HERE",
    url="https://git.ligo.org/ka-lok.lo/skymap-overlap",
    packages=[
        "skymap_overlap",
    ],
    entry_points={
        "console_scripts": [
            "compute_overlap=skymap_overlap.compute_overlap:main",
            "download_skymap=skymap_overlap.download_skymap:main",
            "compute_overlap_from_skymaps_pipe=skymap_overlap.compute_overlap_from_skymaps_pipe:main",
        ]
    },
    install_requires=[
        "healpy",
        "ligo.skymap",
        "pycondor",
    ],
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX",
    ],
    python_requires='>=3.6',
)
