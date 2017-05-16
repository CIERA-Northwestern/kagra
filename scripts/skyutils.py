import copy
import itertools

from scipy.optimize import brentq
import healpy
import numpy as np
from matplotlib import pyplot as plt

from lal import LIGOTimeGPS, GreenwichMeanSiderealTime, ComputeDetAMResponse
import lal
detectors = dict([(d.frDetector.prefix, d) for d in lal.CachedDetectors])

import lalsimulation
"""
_ref_h, _ = lalsimulation.SimInspiralFD(
            1.4 * lal.MSUN_SI,  1.4 * lal.MSUN_SI,
            0., 0., 0.,
            0., 0., 0.,
            100e6 * lal.PC_SI, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.125, 10.0, 2048., None,
            lalsimulation.SimInspiralGetApproximantFromString("IMRPhenomPv2"))
"""
_ref_h_bns, _ = lalsimulation.SimInspiralFD(0.0, 0.125,
            1.4 * lal.MSUN_SI,  1.4 * lal.MSUN_SI,
            0., 0., 0.,
            0., 0., 0.,
            10., 2048., 10.,
            100e6 * lal.PC_SI, 0.0, 0.0,
            0.0, 0.0, None, None, -1, -1,
            lalsimulation.SimInspiralGetApproximantFromString("IMRPhenomPv2"))

_ref_h_nsbh, _ = lalsimulation.SimInspiralFD(0.0, 0.125,
            1.4 * lal.MSUN_SI,  10. * lal.MSUN_SI,
            0., 0., 0.,
            0., 0., 0.,
            10., 2048., 10.,
            100e6 * lal.PC_SI, 0.0, 0.0,
            0.0, 0.0, None, None, -1, -1,
            lalsimulation.SimInspiralGetApproximantFromString("IMRPhenomPv2"))

_ref_h_nsbh_ms, _ = lalsimulation.SimInspiralFD(0.0, 0.125,
            1.4 * lal.MSUN_SI,  10. * lal.MSUN_SI,
            0., 0., 0.,
            0., 0., 0.9,
            10., 2048., 10.,
            100e6 * lal.PC_SI, 0.0, 0.0,
            0.0, 0.0, None, None, -1, -1,
            lalsimulation.SimInspiralGetApproximantFromString("IMRPhenomPv2"))


_default_psds = {
    "H1": lalsimulation.SimNoisePSDaLIGOZeroDetHighPower,
    "I1": lalsimulation.SimNoisePSDaLIGOZeroDetHighPower,
    "K1": lalsimulation.SimNoisePSDaLIGOZeroDetHighPower,
    "L1": lalsimulation.SimNoisePSDaLIGOZeroDetHighPower,
    "V1": lalsimulation.SimNoisePSDAdvVirgo
}

def effective_bandwidth(psd=lalsimulation.SimNoisePSDaLIGOZeroDetHighPower, h=None, get_timing=True):
    """
    Fairhurst 2009, eqn 23.
    """
    flow = int(10 / 0.125)
    freq = np.linspace(0, 2048, int(2048/0.125) + 1)
    freq = freq[flow:]

    if h is None:
        h = _ref_h_bns.data.data[flow:]
    else:
        h = h.data.data[flow:]
    psd = map(psd, freq)
    psd = np.asarray(psd)

    ip = np.real(h.conj() * h / psd) * 0.125
    snr = 2 * np.sqrt(np.sum(ip))

    mean = 4 * np.sum(ip * freq) / snr ** 2
    var = 4 * np.sum(ip * freq**2) / snr ** 2

    sig_f = np.sqrt(var - mean**2)
    if get_timing:
        return snr, (1.0 / 2 / np.pi / sig_f / snr)
    else:
        return sig_f

# Some reference numbers from Fairhurst 2009
# SNR of 10, up to ISCO
def get_timing_dict(psd_dict={}, snr=10, sig=_ref_h_bns):
    _timing = {}
    for d, psd in psd_dict.iteritems():
        if psd is None:
            psd = _default_psds[d]
        _snr, delta_t = effective_bandwidth(psd, sig)
        delta_t *= _snr / snr
        _timing[d] = delta_t

    return _timing

def scale_timing_dict(timing, snr=10, dets=None):
    """
    Get a dictionary with the \sigma_t values for each detector specified by 'dets'. If 'dets' is not specified, all five HKILV will be retrieved. The timing errors are scaled by signal to noise ratio, and pinned to an SNR of 10.
    """
    new_dict = copy.copy(timing)
    for d in copy.copy(new_dict.keys()):
        if dets is not None and d not in dets:
            del new_dict[d]
            continue
        # Timing uncertainty (uncorrelated) scales like 1/snr
        new_dict[d] *= 10. / snr
    return new_dict

def iter_dets(dets, n=4, unique=True):
    #examined = set()
    for comb in itertools.product(*([dets] * n)):
        net = "".join(sorted(set([d[0] for d in comb])))
        if unique and len(set(comb)) != n - 1:
            continue
        if (comb[0] == comb[1]) or (comb[2] == comb[3]):
            continue
        if comb[1] != comb[2]:
            continue
        # FIXME: It seems that the redundancy is actually required
        #if net in examined:
            #continue
        #examined.add(net)
        yield comb

def baseline(d1, d2):
    r1 = detectors[d1].location
    r2 = detectors[d2].location
    return r1 - r2

def source_vec_from_pos(ra, dec):
    """
    RA, DEC -> normalized position vector
    TODO: This probably should be negative.
    """
    return np.asarray((np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)))

def source_pos_from_vec(x, y, z):
    """
    normalized position vector -> RA, DEC
    """
    r = np.sqrt(x**2 + y**2 + z**2)
    th = 0 if r == 0 else np.arccos(z / r)
    phi = np.arctan2(y, x)
    return np.asarray((th, phi))

def root_plane_intercept(plane_vec, r=lal.REARTH_SI):
    # Generate a set of "right ascensions"
    # FIXME: This will be a simgular problem if the plane is parallel with a
    # meridien of longitude, and will be bad for planes which are close to this
    phis = np.linspace(0, np.pi * 2, 1000)
    thetas = []
    for cphi, sphi in zip(np.cos(phis), np.sin(phis)):
        # This is the function we're minimizing
        def planar_distance(theta):
            rad_vec = (r*np.sin(theta)*cphi, r*np.sin(theta)*sphi, r*np.cos(theta), 1)
            return np.dot(plane_vec, rad_vec)
        theta_max = brentq(planar_distance, 0, np.pi)
        thetas.append(theta_max)
    return thetas, phis

def detector_plane(network):
    assert len(network) == 3
    r_1 = baseline(network[0], network[1])
    r_2 = baseline(network[1], network[2])
    det_plane = np.cross(r_1, r_2)
    det_plane /= np.sqrt(np.dot(det_plane, det_plane))

    p0 = detectors[network[0]].location
    itrcp = -np.dot(det_plane, p0)

    return root_plane_intercept(np.append(det_plane, itrcp))

def seconds_since_utc_midnight(gpstime):
    return (gpstime / 3600. % 24.) * 60 * 60

_default_gps_time = 1e9 - seconds_since_utc_midnight(1e9)
def solid_angle_error(source_pos, network_timing, gps_time=_default_gps_time):
    """
    Schutz, 2011, eqn. 31 --- quoted from earlier Wen and Chen.
    """
    source_pos = np.asarray(source_pos)
    del_omega_inv_sq = np.zeros(source_pos.shape[-2:])
    # translate sources to utc midnight
    source_pos[0] -= seconds_since_utc_midnight(gps_time) / 3600. * np.pi / 12
    source_vec = source_vec_from_pos(*source_pos)
    for det_comb in iter_dets(network_timing.keys()):
        timing = map(network_timing.__getitem__, det_comb)
        # NOTE: We reverse the components because of the indexing, but doesn't
        # really matter once we take the magnitude
        r_1 = baseline(det_comb[1], det_comb[0])
        r_2 = baseline(det_comb[3], det_comb[2])
        det_plane = np.cross(r_1, r_2)
        # FIXME: Bad, don't know why
        det_plane[-1] *= -1
        del_omega_inv_sq += np.prod(timing)**-2 * np.dot(det_plane, np.swapaxes(source_vec, 1, 0))**2
    timing_sum = np.sum(1. / np.power(network_timing.values(), 2))

    return 4 * np.sqrt(2) * np.pi * lal.C_SI**2 * timing_sum / np.sqrt(del_omega_inv_sq)

def delay_to_plane(source_pos, network_timing):
    """
    Schutz, 2011, eqn. 31 --- quoted from earlier Wen and Chen.
    """
    del_omega_inv_sq = np.zeros(np.shape(source_pos)[-2:])
    source_vec = source_vec_from_pos(*source_pos)
    for det_comb in iter_dets(network_timing.keys()):
        # NOTE: We reverse the components because of the indexing, but doesn't
        # really matter once we take the magnitude
        r_1 = baseline(det_comb[1], det_comb[0])
        r_2 = baseline(det_comb[3], det_comb[2])
        if np.all(np.cross(r_1, r_2) == 0):
            continue
        r_1 /= np.sqrt(np.dot(r_1, r_1))
        r_2 /= np.sqrt(np.dot(r_2, r_2))
        det_plane = np.cross(r_1, r_2)
        det_plane /= np.sqrt(np.dot(det_plane, det_plane))
        # TODO: Bad, don't know why
        det_plane[-1] *= -1
        ip = []
        _tmp = np.asarray(source_pos).reshape(2, -1)
        for i, vec in enumerate(source_vec.reshape(3, -1).T):
            ip.append(np.dot(det_plane, vec))
        del_omega_inv_sq += np.dot(det_plane, np.swapaxes(source_vec, 1, 0))
        break

    return del_omega_inv_sq

def create_table(fitsfile, verbose=False):
    skymap = healpy.read_map(fitsfile, verbose=verbose)
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
    ra_int = np.linspace(-12, 12., 2*npts, endpoint=False)
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

    # form RA and dec into coordinates expected by basemap
    ra_int *= 12 / np.pi
    ra_int -= 12
    dec_int -= np.pi/2
    dec_int *= -1
    dec_int *= 180 / np.pi

    return ra_int, dec_int, prob_int

def _sph_grid(npts):
    ra_grid = np.linspace(-np.pi, np.pi, 2*npts)
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

