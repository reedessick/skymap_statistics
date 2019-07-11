description = """compute quantities about skymaps based on known detector locations and how triangulation works"""
author = "Reed Essick (reed.essick@ligo.org)"

#=================================================

import numpy as np
import healpy as hp

import detector_cache 

from lal.lal import GreenwichMeanSiderealTime as GMST

#=================================================

rad2deg = 180/np.pi
deg2rad = 1./rad2deg

twopi = 2*np.pi

c = 299792458.0 #m/s

#=================================================

detectors = detector_cache.detectors

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

def rotateRAC2C( ra, gps1, gps2, noWRAP=False ):
    """
    rotates the RA according to the change in gps

    takes ra at gps1 and rotates it so that the earth-fixed coordinates are invarient but the time has changed to gps2
    """
    gmst2 = GMST( gps2 )
    gmst1 = GMST( gps1 )

    if noWRAP:
        return ra - (gmst1-gmst2)%twopi
    else:
        return (ra - gmst1 + gmst2)%(twopi)

def rotateRAC2E( ra, gps, noWRAP=False ):
    """
    rotates ra -> earth fixed coords
    """
    gmst = GMST( gps )
    if noWRAP:
        return ra - gmst%(twopi)
    else:
        return (ra - gmst)%(twopi)

def rotateRAE2C( phi, gps, noWRAP=False ):
    """
    rotates earth fixed coords -> ra
    """
    gmst = GMST( gps )
    if noWRAP:
        return phi + gmst%(twopi)
    else:
        return (phi + gmst)%(twopi)

def rotateMap( posterior, dphi ):
    """
    rotates phi -> phi+dphi
    """
    npix = len(posterior)
    nside = hp.npix2nside( npix )
    theta, phi = hp.pix2ang( nside, np.arange(npix) )
    phi += dphi

    new_pix = hp.ang2pix( nside, theta, phi )

    return posterior[new_pix]

def rotateMapC2C( posterior, old_gps, new_gps ):
    """
    returns a rotated map that keeps the posterior in the same relative position to the detectors at the new_gps time 
    as it was at the old_gps time.
    """
    npix = len(posterior)
    nside = hp.npix2nside( npix )

    theta, new_phi = hp.pix2ang( nside, np.arange( npix ) )
    phi = rotateRAC2C( new_phi, new_gps, old_gps ) ### rotate the RA according to times
                                                   ### need to map the RAs at the new gps time into the RAs at the old gps time

    new_pix = hp.ang2pix( nside, theta, phi )

    return posterior[new_pix]

def rotateMapC2E( posterior, gps ):
    npix = len(posterior)
    nside = hp.npix2nside( npix )

    theta, phi = hp.pix2ang( nside, np.arange( npix ) )
    ra = rotateRAE2C( phi, gps ) ### rotate phi to get ra -> needed to determine original indexing

    new_pix = hp.ang2pix( nside, theta, ra )

    return posterior[new_pix]

def rotateMapE2C( posterior, gps ):
    npix = len(posterior)
    nside = hp.npix2nside( npix )

    theta, ra = hp.pix2ang( nside, np.arange( npix ) )
    phi = rotateRAC2E( ra, gps ) ### rotate the RA to get phi -> needed to determine original indexing

    new_pix = hp.ang2pix( nside, theta, phi )

    return posterior[new_pix]


#=================================================

def antipode( x, y, coord="C", degrees=False ):
    """
    if coord=C : x->ra, y->dec
    if coord=E : x->phi, y->theta
    returns the anitpode position
    """
    if coord=="C":
        Y, X = hp.vec2ang( -hp.ang2vec(0.5*np.pi-y,x) )
        Y = 0.5*np.pi - Y
    elif coord=="E":
        Y, X = hp.vec2ang( -hp.ang2vec(y,x) )
    return X, Y

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
    dr = detectors[ifo2].dr - detectors[ifo1].dr
    if normed:
        dr /= np.sum(dr**2)**0.5
    return dr

#=================================================

def time_delay( theta, phi, ifo1, ifo2, coord="E", tgeocent=None, degrees=False ):
    """
    returns the time delay between two ifo's (in seconds)

    interprets direction as (theta, phi) or (dec, ra) depending on \"coord\"
    """
    if ifo1 not in detectors.keys():
        raise ValueError("ifo1=%s not understood"%ifo1)
    if ifo2 not in detectors.keys():
        raise ValueError("ifo2=%s not understood"%ifo2)

    dr = detectors[ifo1].dr - detectors[ifo2].dr

    if degrees:
        theta *= deg2rad
        phi *= deg2rad

    if coord=="C":
        theta, phi = __celest2earth( theta, phi, tgeocent )

    theta = theta%np.pi
    phi = phi%(2*np.pi)

    if isinstance(theta, float):
        return np.sum( hp.ang2vec( theta, phi)*dr )
    else:
        return np.sum( hp.ang2vec( theta, phi)*dr, axis=1 )

def time_delay_locus( dt, ifo1, ifo2, coord="E", tgeocent=None, Nsamp=1001, degrees=False ):
    """
    returns a locus of points with the time delay between ifo's equal to dt
    """
    if ifo1 not in detectors.keys():
        raise ValueError("ifo1=%s not understood"%ifo1)
    if ifo2 not in detectors.keys():
        raise ValueError("ifo2=%s not understood"%ifo2)

    ### compute associated polar angle in line-of-sight frame 
    dr = detectors[ifo1].dr - detectors[ifo2].dr
    cosTheta = dt / np.sum(dr*dr)**0.5

    ### in line-of-sight frame, define curves
    theta = np.ones((Nsamp,), dtype=float)*np.arccos(cosTheta)
    phi = np.linspace(0, 2*np.pi, Nsamp)

    ### rotate from line-of-sight frame to Earth-Fixed frame
    thetaPole, phiPole = line_of_sight( ifo2, ifo1, coord="E", tgeocent=None, degrees=False )
    theta, phi = rotate2pole( theta, phi, -thetaPole, phiPole, degrees=False )

    ### rotate if needed and convert to degrees if needed
    if coord=="C":
        dec, ra = __earth2celest( theta, phi, tgeocent )
        
        if degrees:
            dec *= rad2deg
            phi *= rad2deg
        return dec, ra

    elif coord=="E":
        if degrees:
            theta *= rad2deg
            phi *= rad2deg
        return theta, phi

#=================================================

def overhead( ifo, coord="E", tgeocent=None, degrees=False ):
    """
    get the overhead direction of an ifo. 

    overhead is defined as (nx x ny) = nz
    """
    if ifo not in detectors.keys():
        raise ValueError("ifo=%s not understood"%ifo )

    t, p = __overhead_Evec( ifo )

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
    return detectors[ifo].zenith

#=================================================

def rotate2pole( theta, phi, thetaPole, phiPole, degrees=False ):
    """
    rotates the set of coordinates (theta, phi) into a new corrdinate system
    (thetaPole, phiPole) define the new coordinate pole, with phi=phiPole corresonding to the new zero of the azimuth.
    accomplished through successive applications of healpy.Rotator objects
    """
    newt, newp = hp.Rotator( deg=degrees, rot=[phiPole, -thetaPole, 0], eulertype="ZYZ" )( theta, phi )
    newt, newp = hp.Rotator( deg=degrees, rot=[-phiPole], eulertype="Z" )( newt, newp )

    return newt, newp

#=================================================

def mutualinformation( count, bins=None ):
    """
    returns the KL distance between the joint distribution and the product of the marginals
        sum p(t,p) * log( p(t,p)/p(t)*p(p) )
    if bins!=None:
        theta_bins, phi_bins = bins
        include the jacobian for spherical coordinates in the summation (approximate an integral instead of a sum)

    assumes count is indexed by phi, then theta
    -> p(theta) = np.sum( count, axis=0 )
    """
    if bins:
        theta, phi = bins
        weights = np.outer( np.cos(theta)[:-1]-np.cos(theta)[1:] , phi[1:]-phi[:-1] ) ### weight by the area contained in that bin
    else:
        weights = np.ones_like( count ) ### weight all bins equally

    ### ensure normalization
    count /= np.sum( weights*count )

    ### compute marginals
    marg0 = np.sum( weights*count, axis=1 ) / np.sum( weights, axis=1 )
    marg0 /= np.sum( marg0 )
    marg1 = np.sum( weights*count, axis=0 ) / np.sum( weights, axis=0 )
    marg1 /= np.sum( marg1 )

    margs = np.outer( marg0, marg1 )
    margs /= np.sum(margs)

    truth = weights*count > 0

    mi =  np.sum( (weights*count)[truth] * np.log( (weights*count)[truth] / margs[truth] ) )

#    truth = marg0 > 0
#    ent0 = np.sum( marg0[truth] * np.log( marg0[truth] ) )

#    truth = marg1 > 0
#    ent1 = np.sum( marg1[truth] * np.log( marg1[truth] ) )

    truth = count > 0
    entj = np.sum( (weights*count)[truth] * np.log( (weights*count)[truth] ) )

    return mi, -entj

def compute_mi( theta, phi, Nbins, weights=None ):
    theta_bins = np.linspace(0, np.pi, Nbins+1)
    phi_bins = np.linspace(-np.pi, np.pi, Nbins+1)

    count = np.histogram2d( phi, theta, bins=(phi_bins, theta_bins), weights=weights )[0]

    return mutualinformation( count, bins=None )

