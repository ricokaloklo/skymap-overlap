import setuptools
from pathlib import Path

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
    author_email="ka-lok.lo@ligo.org",
    description="Compute overlap between two skymaps",
    long_description=Path("README.md").read_text(encoding="utf-8"),
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
    python_requires='>=3.9',
)
