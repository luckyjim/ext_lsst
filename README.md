# ext_lsst
LSST for EUCLID

# Initialisation package


Go to root package

```console
cd /path/to/EXT_LSST_Testing
```
and source init file

```console
. ./init_stage1.bash
```
# Package presentation
## Description of "at_ccin2p3" directory

Perform the LSST stage 1 in 2 step "LSST env" and "Euclid env" with scheduler of the CCIN2P3. See [readme.md](at_ccin2p3/readme.md) for more informations

## Description of "convert_tools" directory
Mainly:
* convert WCS SIP (used by LSST) to WCS TPV (used by Euclid)
* convert Euclid reference stars catalog (True Univers for simulation, GAIA for true data) to LSST stack format (texte format )

## Description of "EXT_PF1_dm4lsst" directory

This directory has the "official" Euclid/Elements structure to provide xml file (Data model) associated to FITS product provide by LSST stage 1. See [readme.md](EXT_PF1_dm4lsst/readme.md) for more informations
