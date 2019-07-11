description = " a holder for default/known detectors "
author      = "reed.essick@ligo.org"

#=================================================

import copy
import antenna

from pkg_resources import resource_filename

import numpy as np

#=================================================

class PSD(object):
    """
    an object that holds onto power-spectral densities with associated frequency samples
    we define a scipy.interpolate.interp1d object for convenience
    """

    def __init__(self, freqs, psd, kind="linear"):
        len_freqs = len(freqs)
        self.n_freqs = len_freqs
        if len(psd) != len_freqs:
            raise ValueError, "freqs and ps must have the same length"
        if not len_freqs:
            raise ValueError, "freqs and psd must have at least 1 entries"
        elif len_freqs == 1:
            freqs = np.array(2*list(freqs))
            psd = np.array(2*list(psd))
        self.freqs = freqs
        self.psd = psd

    def check(self):
        return len(self.freqs) == len(self.psd)

    def update(self, psd, freqs=None):
        if freqs!=None:
            if len(freqs)!=len(psd):
                raise ValueError, "len(freqs) != len(psd)"
            self.freqs = freqs[:]
            self.psd = psd[:]
        else:
            self.psd=psd[:]

    def get_psd(self):
        return self.psd

    def get_freqs(self):
        return self.freqs

    def interpolate(self, freqs):
        return np.interp(freqs, self.freqs, self.psd)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        min_psd = np.min(self.psd)
        d=int(np.log10(min_psd))-1
        return """PSD object
        min{freqs}=%.5f
        max{freqs}=%.5f
        No. freqs =%d
        min{psd}=%.5fe%d  at freqs=%.5f"""%(np.min(self.freqs), np.max(self.freqs), len(self.freqs), min_psd*10**(-d), d, self.freqs[min_psd==self.psd][0])

class Detector(object):
    """
    an object representing a gravitational wave detector. methods are meant to be convenient wrappers for more general operations. 
    """

    def __init__(self, name, dr, nx, ny, psd):
        """
        name = None  # detector's name (eg: H1)
        dr = np.zeros((3,)) # r_detector - r_geocent
        nx = np.zeros((3,)) # direction of the x-arm
        ny = np.zeros((3,)) # direction of the y-arm
        psd = None   # the psd for network (should be power, not amplitude)
        """
        self.name = name
        if not isinstance(dr, np.ndarray):
            dr = np.array(dr)
        self.dr = dr
        if not isinstance(nx, np.ndarray):
            nx = np.array(nx)
        self.nx = nx
        if not isinstance(ny, np.ndarray):
            ny = np.array(ny)
        self.ny = ny

        self.psd = psd

    def __str__(self):
        return "Detector : %s"%self.name

    def __repr__(self):
        return self.__str__()

    def set_psd(self, psd, freqs=None):
        self.psd.update(psd, freqs=freqs)

    def get_psd(self):
        return self.psd

    def antenna_patterns(self, theta, phi, psi, freqs=None, dt=None):
        """ returns the antenna patterns for this detector. If psi is not supplied, returns antenna patterns that diagonalize A_{ij} """
        if dt != None:
            return antenna.antenna_patterns(theta, phi, psi, self.nx, self.ny, freqs=freqs, dt=dt)
        else:
            return antenna.antenna_patterns(theta, phi, psi, self.nx, self.ny, freqs=freqs, dr=self.dr)

    def project(self, theta, phi, psi, hp, hx, freqs=None):
        Fp, Fx = self.antenna_patterns(theta, phi, psi, freq=freqs)
        return Fp*hp + Fx*hx

    def snr(self, data, freqs=None):
        """ 
        returns the SNR for data using the PSD stored within the object 
        if freqs==None: assumes data corresponds to self.psd.freqs
        """
        if freqs==None:
            freqs = self.get_psd().get_freqs()
        if len(data) != len(freqs):
            raise ValueError, "len(data) != len(freqs)"

        return ( 4*np.sum((data.real**2+data.imag**2) / self.get_psd().interpolate(freqs))*(freqs[1]-freqs[0]) )**0.5 ### return SNR

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return """Detector object
        name : %s
        dr = %.5f , %.5f , %.5f 
        nx = %.5f , %.5f , %.5f
        ny = %.5f , %.5f , %.5f
        PSD : %s"""%(self.name, self.dr[0], self.dr[1], self.dr[2], self.nx[0], self.nx[1], self.nx[2], self.ny[0], self.ny[1], self.ny[2], str(self.psd))

class Network(object):
    """an object respresenting a network of detectors
    """
    def __init__(self, detectors=[]):
        self._detectors = dict((det.name, det) for det in detectors)

    def __len__(self):
        return len(self._detectors)

    @property
    def _instr(self):
        return list(self._detectors.keys())

    def add(self, det):
        assert det.name not in self._detectors, 'detector %s already included in the network!'%det.name
        self._detectors[det.name] = det

    def remove(self, det):
        if not isinstance(det, str):
            det = det.name
        assert det in self._detectors, 'detector %s not in the network!'%det
        self._detectors.pop(det)

    def snr(self, *args, **kwargs):
        """a wrapper that will attempt to automatically detect the function signature and delegate as needed
        """
        if len(args)==5: ### interpret as (theta, phi, psi, hp, hx)
            theta, phi, psi, hp, hx = args
            freqs = kwargs.get('freqs', None)

        elif len(args)==1: ### interpret as (event,)
            event = args[0]
            theta = event.theta
            phi = event.phi
            psi = event.polarization
            hp, hx, freqs = event.waveform(
                flow=kwargs.get('flow', 10.),
                fhigh=kwargs.get('fhigh', 2048.),
                deltaf=kwargs.get('deltaf', 0.125),
            )

        else:
            raise NotImplementedError('function signature not recognized. Please specify either (theta, phi, psi, hp, hx) or (event,)')

        return dict((det.name, det.snr(det.project(theta, phi, psi, hp, hx, freqs=freqs))) for det in self._detectors.values())

#=================================================
# known PSDs
#=================================================
psds = {}

default_psd = PSD(np.array([0]), np.array([1]), kind="linear")
psds["default"] = default_psd

### design psd's
aligo_design_psd_file = resource_filename(__name__, 'PSDs/aLIGO_design.txt')
aligo_design_psd_dat = np.genfromtxt(aligo_design_psd_file)
aligo_design_psd = PSD(aligo_design_psd_dat[:,0], aligo_design_psd_dat[:,1]**2, kind="linear")
psds["aligo_design"] = aligo_design_psd

avirgo_design_psd_file = resource_filename(__name__, 'PSDs/aVirgo_design.txt')
avirgo_design_psd_dat = np.genfromtxt(avirgo_design_psd_file)
avirgo_design_psd = PSD(avirgo_design_psd_dat[:,0], avirgo_design_psd_dat[:,1]**2, kind="linear")
psds["avirgo_design"] = avirgo_design_psd

#=================================================
# known detectors
#=================================================

c = 299792458.0 #m/s

### Detector locations and orientations taken from Anderson, et all PhysRevD 63(04) 2003
detectors = {}

__H_dr__ = np.array((-2.161415, -3.834695, +4.600350))*1e6/c # sec
__H_nx__ = np.array((-0.2239, +0.7998, +0.5569))
__H_ny__ = np.array((-0.9140, +0.0261, -0.4049))
detectors["H"] = Detector("H", __H_dr__, __H_nx__, __H_ny__, copy.deepcopy(aligo_design_psd))

__L_dr__ = np.array((-0.074276, -5.496284, +3.224257))*1e6/c # sec
__L_nx__ = np.array((-0.9546, -0.1416, -0.2622))
__L_ny__ = np.array((+0.2977, -0.4879, -0.8205))
detectors["L"] = Detector("L", __L_dr__, __L_nx__, __L_ny__, copy.deepcopy(aligo_design_psd))

__V_dr__ = np.array((+4.546374, +0.842990, +4.378577))*1e6/c # sec
__V_nx__ = np.array((-0.7005, +0.2085, +0.6826))
__V_ny__ = np.array((-0.0538, -0.9691, +0.2408))
detectors["V"] = Detector("V", __V_dr__, __V_nx__, __V_ny__, copy.deepcopy(avirgo_design_psd))

