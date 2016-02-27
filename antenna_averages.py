#!/usr/bin/python
usage = "anteanna_averages.py [--options] FITS [gps]"
description = "a simple script that computes averages of a few antenna pattern things. Relies on bayesburst code"
author = "reed.essick@ligo.org"

import numpy as np
import healpy as hp

import detector_cache ### relies on BayesBurst code

from optparse import OptionParser

#-------------------------------------------------

parser = OptionParser(usage=usage, description=description)

parser.add_option("-v", "--verbose", default=False, action="store_true")

parser.add_option("-o", "--observatory", default=[], type="string", action="append")
parser.add_option("-c", "--coord", default="C", type="string", help="either E or C")

parser.add_option("", "--gps", default=None, type="float", help="must be supplied if --coord=C")

opts, args = parser.parse_args()

if opts.coord=="C" and opts.gps==None:
    opts.gps = float(raw_input("gps = "))

if len(args)!=1:
    raise ValueError("Please supply exactly one input argument\n%s"%(usage))
fitsname = args[0]

#-------------------------------------------------

if opts.verbose:
    print "reading map from %s"%fitsname

### load map
map = hp.read_map( fitsname )
mapind = map.argmax()

npix = len(map)
nside = hp.npix2nside( npix )

theta, phi = hp.pix2ang( nside, np.arange(npix) )
psi = 0.0 ### may want to pick a better choice than this...

if opts.coord=="C": ### rotate ra->phi
    from lal.lal import GreenwichMeanSiderealTime as GMST
    gmst = GMST( opts.gps )
    phi = (phi-gmst)%(2*np.pi)

### find MAP coordinate
mapind = map.argmax()
thetaMAP = theta[mapind]
phiMAP = theta[mapind]

### iterate through observatories
for ifo in opts.observatory:
    if not detector_cache.detectors.has_key(ifo):
        raise ValueError("--observatory=%s not understood"%(ifo))    

    print ifo

    detector = detector_cache.detectors[ifo]

    Fp, Fx = detector.antenna_patterns( theta, phi, psi )
    print "    MAP:\n        F+^2 + Fx^2 = %.9f\n        F+ = %.9f\n        Fx = %.9f"%(Fp[mapind]**2 + Fx[mapind]**2, Fp[mapind], Fx[mapind])
    print "    ave:\n        F+^2 + Fx^2 = %.9f\n        F+ = %.9f\n        Fx = %.9f"%(np.sum( map*(Fp**2 + Fx**2) ), np.sum( map*Fp ), np.sum( map*Fx ))
