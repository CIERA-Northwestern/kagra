import healpy
import numpy as np
from matplotlib import pyplot as plt
from lal import LIGOTimeGPS, GreenwichMeanSiderealTime, ComputeDetAMResponse

import lal
detectors = lal.CachedDetectors
detectors = dict([(d.frDetector.prefix, d) for d in detectors])

def create_table(fitsfile):
    skymap = healpy.read_map(fitsfile)
    dat = np.zeros((4, len(skymap)))
    ns = healpy.npix2nside(len(skymap))
    dec, ra = healpy.pix2ang(ns, np.asarray(range(len(skymap))))
    dat[0,:] = ra * 12 / np.pi
    dat[1,:] = dec * 180 / np.pi
    dat[2,:] = skymap
    dat = dat[:,np.argsort(dat[2,:])[::-1]]
    dat[3,:] = np.cumsum(dat[2,:])
    dat = np.rec.fromarrays(dat, names=("ra", "dec", "prob", "cumul"))
    return dat, skymap

def interpolate_map(ra, dec, prob, npts=100):
    ra_int = np.linspace(0, 24., 2*npts, endpoint=False)
    dec_int = np.linspace(-90, 90, npts, endpoint=False)

    from matplotlib.mlab import griddata
    prob_int = griddata(ra, dec, prob, ra_int, dec_int, interp='linear')
    ra_int, dec_int = np.meshgrid(ra_int, dec_int)
    return ra_int, dec_int, prob_int

def interpolate_healpix_map(ra, dec, prob, npts=100):
    ra_int = np.linspace(0, 2*np.pi, 2*npts)
    dec_int = np.linspace(0, np.pi, npts)

    ra_int, dec_int = np.meshgrid(ra_int, dec_int)
    prob_int = healpy.pixelfunc.get_interp_val(prob, dec_int, ra_int)

    ra_int *= 12 / np.pi
    dec_int -= np.pi/2
    dec_int *= 180 / np.pi

    return ra_int, dec_int, prob_int

def _sph_grid(npts):
    ra_grid = np.linspace(0, 2*np.pi, 2*npts)
    dec_grid = np.linspace(-np.pi/2, np.pi/2, npts)
    return np.meshgrid(ra_grid, dec_grid)

def net_antenna_pattern(gpstime, network, psi=0, npts=100, norm=False):

    # FIXME: need to check on this
    gps = LIGOTimeGPS(gpstime)
    gmst_rad = GreenwichMeanSiderealTime(gps)

    ra_grid, dec_grid = _sph_grid(npts)

    net_pat = np.zeros(ra_grid.shape[0]*ra_grid.shape[1])

    net_align = np.zeros(ra_grid.shape[0]*ra_grid.shape[1])
    net_dpf = np.zeros(ra_grid.shape[0]*ra_grid.shape[1])

    psi_rad = 0
    i = 0
    for ra_rad, de_rad in zip(ra_grid.flat, dec_grid.flat):
        fp = [ComputeDetAMResponse(detectors[ifo].response, ra_rad, de_rad, psi_rad, gmst_rad)[0] for ifo in network]
        fx = [ComputeDetAMResponse(detectors[ifo].response, ra_rad, de_rad, psi_rad, gmst_rad)[1] for ifo in network]
        fp = np.asarray(fp)
        fx = np.asarray(fx)
        fp2, fx2 = np.dot(fp, fp), np.dot(fx, fx)

        net_dpf[i] = psi_dpf = 0.5 * np.arctan2(2 * np.dot(fp, fx), (fp2 - fx2))

        fp, fx = fp * np.cos(psi_dpf) + fx * np.sin(psi_dpf), \
                -fp * np.sin(psi_dpf) + fx * np.cos(psi_dpf)
        fp2, fx2 = np.dot(fp, fp), np.dot(fx, fx)
        net_pat[i] = np.sqrt(fp2 + fx2)
        net_align[i] = np.sqrt(fx2 / fp2)
        i += 1

    ra_grid *= 180/np.pi
    dec_grid *= 180/np.pi

    net_pat = net_pat.reshape(ra_grid.shape)
    net_align = net_align.reshape(ra_grid.shape)
    net_dpf = net_dpf.reshape(ra_grid.shape) / 2 # we multiplied by two above

    if norm:
        net_pat /= np.sqrt(2) #* len(network)

    return ra_grid, dec_grid, net_pat, net_align, net_dpf


def net_antenna_pattern_point(gpstime, network, ra_rad, de_rad, psi=0, npts=100, norm=False):

    """
    Only get the network antenna pattern at a given ra and dec of interest
    """

    # FIXME: need to check on this
    gps = LIGOTimeGPS(gpstime)
    gmst_rad = GreenwichMeanSiderealTime(gps)

    psi_rad = 0
    #i = 0

    fp = [ComputeDetAMResponse(detectors[ifo].response, ra_rad, de_rad, psi_rad, gmst_rad)[0] for ifo in network]
    fx = [ComputeDetAMResponse(detectors[ifo].response, ra_rad, de_rad, psi_rad, gmst_rad)[1] for ifo in network]
    #print network
    #for ifo in network:
    #    print ifo, ComputeDetAMResponse(detectors[ifo].response, ra_rad, de_rad, psi_rad, gmst_rad)[0]

    fp = np.asarray(fp)
    fx = np.asarray(fx)
    fp2, fx2 = np.dot(fp, fp), np.dot(fx, fx)

    net_dpf = psi_dpf = 0.5 * np.arctan2(2 * np.dot(fp, fx), (fp2 - fx2))

    fp, fx = fp * np.cos(psi_dpf) + fx * np.sin(psi_dpf), \
            -fp * np.sin(psi_dpf) + fx * np.cos(psi_dpf)
    fp2, fx2 = np.dot(fp, fp), np.dot(fx, fx)
    net_pat = np.sqrt(fp2 + fx2)
    net_align = np.sqrt(fx2 / fp2)

    #net_pat = net_pat.reshape(ra_grid.shape)
    #net_align = net_align.reshape(ra_grid.shape)
    #net_dpf = net_dpf.reshape(ra_grid.shape) / 2 # we multiplied by two above

    if norm:
        net_pat /= np.sqrt(2) #* len(network)

    return net_pat, net_align, net_dpf





def net_gradient(net_pat, npts=20):
    dra = np.pi/npts
    ddec = np.pi/npts
    grad_net_pat = np.gradient(net_pat, dra, ddec)

    return grad_net_pat

