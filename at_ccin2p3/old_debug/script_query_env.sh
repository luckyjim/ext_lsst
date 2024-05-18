#!/bin/bash

source /cvmfs/sw.lsst.eu/linux-x86_64/lsst_distrib/w_2018_31/loadLSST.bash
setup lsst_distrib

echo "conda is: " $(which conda)

echo "conda info:"
conda info

echo "python is: " $(which python)

echo "Python version:" $(python --version 2>&1)

echo "numpy version: " $(python -c "import numpy as np; print(np.version.version)")

echo "astropy version: " $(python -c "import astropy; print(astropy.__version__)")

echo "ingestSimImages.py is:" $(which ingestSimImages.py)

