#!/usr/bin/python
usage       = "anteanna_averages.py [--options] FITS"
description = "a simple script that computes averages of a few antenna pattern things. Relies on bayesburst code"
author      = "R. Essick (reed.essick@ligo.org)"

#-------------------------------------------------

import numpy as np
import healpy as hp

import detector_cache ### relies on BayesBurst code

from optparse import OptionParser

#-------------------------------------------------

def summarize( post, ifo, coord='C', gps=None, fitsname=None ):
    '''
    computes the summary information about antenna patterns for post
    '''
    npix = len(post)
    nside = hp.npix2nside( npix )

    theta, phi = hp.pix2ang( nside, np.arange(npix) )
    psi = 0.0 ### may want to pick a better choice than this...

    if coord=="C": ### rotate ra->phi
        from lal.lal import GreenwichMeanSiderealTime as GMST
        gmst = GMST( opts.gps )
        phi = (phi-gmst)%(2*np.pi)

    ### find MAP coordinate
    mapind = post.argmax()
    thetaMAP = theta[mapind]
    phiMAP = theta[mapind]

    if not detector_cache.detectors.has_key(ifo):
        raise ValueError("ifo=%s not understood"%(ifo))

    Fp, Fx = detector_cache.detectors[ifo].antenna_patterns( theta, phi, psi )
    if fitsname:
        hp.write_map( fitsname, Fp**2+Fx**2, column_names='Fx^2+F+^2' )

    return (Fpind]**2 + Fx[mapind]**2, p[mapind], Fx[mapind]), (np.sum( map*(Fp**2 + Fx**2) ), np.sum( map*Fp ), np.sum( map*Fx ))

#-------------------------------------------------

if __name__=="__main__":

    parser = OptionParser(usage=usage, description=description)

    parser.add_option("-v", "--verbose", default=False, action="store_true")

    parser.add_option("-o", "--observatory", default=[], type="string", action="append")
    parser.add_option("-c", "--coord", default="C", type="string", help="either E or C")

    parser.add_option("", "--gps", default=None, type="float", help="must be supplied if --coord=C")

    parser.add_option("", "--write-fits", default=False, action='store_true', help='write F+^2+Fx^2 to disk')
    parser.add_option("", "--output-dir", default=".", type='string')
    parser.add_option("", "--tag", default="", type="string")

    opts, args = parser.parse_args()

    if opts.coord=="C" and opts.gps==None:
        opts.gps = float(raw_input("gps = "))

    if len(args)!=1:
        raise ValueError("Please supply exactly one input argument\n%s"%(usage))
    fitsname = args[0]

    if opts.write_fits and (not os.path.exists(opts.output_dir)):
        os.makedirs(opts.output_dir)

    if opts.tag:
        opts.tag = "_"+opts.tag

    #---------------------------------------------

    if opts.verbose:
        print "reading map from %s"%fitsname

    ### load map
    post = hp.read_map( fitsname )

    for ifo in opts.observatory:
        if opts.write_fits:
            new_fitsname = "%s/%s-antennaNorm%s.fits"%(opts.output_dir, ifo, opts.tag)
            if opts.verbose:
                print "writing antenna patterns into : %s"%(new_fitsname)
        else:
            new_fitsname = None

        MAP, AVE = summarize( post, ifo, coord=opts.coord, gps=opts.gps, fitsname=new_fitsname )

        print ifo
        print "    MAP:\n        F+^2 + Fx^2 = %.9f\n        F+ = %.9f\n        Fx = %.9f"%MAP
        print "    ave:\n        F+^2 + Fx^2 = %.9f\n        F+ = %.9f\n        Fx = %.9f"%AVE
