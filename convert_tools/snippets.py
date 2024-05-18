import lsst.daf.persistence


path_in = '/sps/euclid/Users/colley/lsst/butler/stage1_4T2_100009/input/'
path_out = '/sps/euclid/Users/colley/lsst/butler/stage1_4T2_100009/output/'

dataid = {'visit': 100009, 'raft': '1,1', 'sensor': '1,1'}

butler_out = lsst.daf.persistence.Butler(path_out)
calexp = butler_out.get('calexp', dataid)

mask = calexp.getMask()
variance = calexp.getVariance()
wcs = calexp.getWcs()
wcs_dict = wcs.getFitsMetadata().toDict()

# Polygon footprint
width = calexp.getWidth()
height = calexp.getHeight()

polygon = [
    wcs.pixelToSky(0, 0),
    wcs.pixelToSky(width, 0),
    wcs.pixelToSky(width, height),
    wcs.pixelToSky(0, height)
]


background = butler_out.get('calexpBackground', dataid)
background_image = background.getImage()

butler_in = lsst.daf.persistence.Butler(path_in)
rawimage = butler_in.get('eimage', dataid)


# What must be written in the FITS file
# filename: EUC_EXT_DETRENDEDFRAME_LSST-{visitid}-{nn}_{date}.fits



