description = """compute quantities about skymaps based on known detector locations and how triangulation works"""
author = "Reed Essick (reed.essick@ligo.org)"

import numpy as np
import healpy as hp

#=================================================

rad2deg = 180/np.pi

c = 299792458.0 #m/s

#=================================================

detectors = {}

__H_dr__ = np.array((-2.161415, -3.834695, +4.600350))*1e6/c # sec
__H_nx__ = np.array((-0.2239, +0.7998, +0.5569))
__H_ny__ = np.array((-0.9140, +0.0261, -0.4049))
detectors["H"] = {"dr":__H_dr__, "nx":__H_nx__, "ny":__H_ny__}

__L_dr__ = np.array((-0.074276, -5.496284, +3.224257))*1e6/c # sec
__L_nx__ = np.array((-0.9546, -0.1416, -0.2622))
__L_ny__ = np.array((+0.2977, -0.4879, -0.8205))
detectors["L"] = {"dr":__L_dr__, "nx":__L_nx__, "ny":__L_ny__}

__V_dr__ = np.array((+4.546374, +0.842990, +4.378577))*1e6/c # sec
__V_nx__ = np.array((-0.7005, +0.2085, +0.6826))
__V_ny__ = np.array((-0.0538, -0.9691, +0.2408))
detectors["V"] = {"dr":__V_dr__, "nx":__V_nx__, "ny":__V_ny__}

#=================================================

def __earth2celest( t, p, tgeocent ):
    from lal.lal import GreenwichMeanSiderealTime as GMST
    gmst = GMST( tgeocent )
    ra = (p+gmst)%(2*np.pi) ### rotate to get RA
    dc = 0.5*np.pi - t
    return dc, ra

def __celest2earth( dec, ra, tgeocent ):
    from lal.lal import GreenwichMeanSiderealTime as GMST
    gmst = GMST( tgeocent )
    p = (ra-gmst)%(2*np.pi) ### rotate to get RA
    t = 0.5*np.pi - dec
    return t, p

#=================================================

def line_of_sight( ifo1, ifo2, coord="E", tgeocent=None, degrees=False ):
    """
    returns the line-of-sight between two detectors in either Earth-fixed ("E") or Celestial ("C") coordinates.
    If coord=="C": we compute the gmst from the geocenter time supplied. If no time is supplied, an error is raised

    line-of-sight is measured from ifo1 to ifo2

    return direction as (theta, phi) or (dec, ra) depending on \"coord\", and angles are measured in radians or degrees depending on \"degrees\"
    """
    if ifo1 not in detectors.keys():
        raise ValueError("ifo1=%s not understood"%ifo1)
    if ifo2 not in detectors.keys():
        raise ValueError("ifo2=%s not understood"%ifo2)

    dr = __line_of_sight_Evec( ifo1, ifo2, normed=True )
    t, p = hp.vec2ang( dr )

    if isinstance(t, np.ndarray):
        t = t[0]
        p = p[0]

    if coord=="E":
        if degrees:
            t *= rad2deg
            p *= rad2deg
        return t, p

    elif coord=="C":
        if tgeocent==None:
            raise ValueError("please supply tgeocent when coord=C")
        dc, ra = __earth2celest( t, p, tgeocent )
        if degrees:
            dc *= rad2deg
            ra *= rad2deg
        return dc, ra
    else:
        raise ValueError("coord=%s not understood"%coord)

def __line_of_sight_Evec( ifo1, ifo2, normed=True ):
    """
    a helper function to get the line of sight as a cartesian vector in Earth-fixed coordinates

    line-of-sight is measured from ifo1 to ifo2
    """
    dr = detectors[ifo2]['dr'] - detectors[ifo1]['dr']
    if normed:
        dr /= np.sum(dr**2)**0.5
    return dr

#=================================================

def overhead( ifo, coord="E", tgeocent=None, degrees=False ):
    """
    get the overhead direction of an ifo. 

    overhead is defined as (nx x ny) = nz
    """
    if ifo not in detectors.keys():
        raise ValueError("ifo=%s not understood"%ifo )

    nz = __overhead_Evec( ifo )
    t, p = hp.vec2ang( nz )

    if isinstance(t, np.ndarray):
        t = t[0]
        p = p[0]

    if coord == "E":
        if degrees:
            t *= deg2rad
            p *= deg2rad
        return t, p
    elif coord == "C":    
        if tgeocent==None:
            raise ValueError("please supply tgeocent when coord=C")
        dc, ra = __earth2celest( t, p, tgeocent )
        if degrees:
            dc *= rad2deg
            ra *= rad2deg
        return dc, ra

def __overhead_Evec( ifo ):
    """
    a helper function to get the overhead direction
    """
    nx = detectors[ifo]['nx']
    ny = detectors[ifo]['ny']
    return np.array( [ nx[1]*ny[2] - nx[2]*ny[1] , -nx[0]*ny[2] + nx[2]*ny[0] , nx[0]*ny[1] - nx[1]*ny[0] ] ) ### take cross product by hand...

#=================================================

def rotate2pole( theta, phi, thetaPole, phiPole, degrees=False ):
    """
    rotates the set of coordinates (theta, phi) into a new corrdinate system
    (thetaPole, phiPole) define the new coordinate pole, with phi=phiPole corresonding to the new zero of the azimuth.
    accomplished through successive applications of healpy.Rotator objects
    """
    newt, newp = hp.Rotator( deg=degrees, rot=[phiPole, -thetaPole, 0], eulertype="ZYZ" )( theta, phi )
    return newt, newp

#=================================================

