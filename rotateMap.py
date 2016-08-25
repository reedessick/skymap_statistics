#!/usr/bin/python
usage = "rotateMap.py [--options] gps input.fits.gz output.fits.gz"
description = "reads in the map and rotates it from geographic to equatorial coordinates or vice versa"
author = "reed.essick@ligo.org"

#-------------------------------------------------

import triangulate
import healpy as hp
import numpy as np

from optparse import OptionParser

#-------------------------------------------------

parser = OptionParser(usage=usage, description=description)

parser.add_option('-v', '--verbose', default=False, action='store_true')

parser.add_option('-S', '--source-coord', default='C', type='string', help='coordinates of the source fits file (either "C" or "E"). DEFAULT="C"')
parser.add_option('-T', '--target-coord', default='E', type='string', help='coordinates of the target fits file (either "C" or "E"). DEFAULT="E"')

opts, args = parser.parse_args()

if len(args)!=3:
    raise ValueError('please supply exactly 3 input arguments\n%s'%usage)
gps = float(args[0])
inFITS  = args[1]
outFITS = args[2]

if opts.source_coord not in ['C', 'E']:
    raise ValueError('--source-coord must be either "C" or "E"')
if opts.target_coord not in ['C', 'E']:
    raise ValueError('--target-ccord must be either "C" or "E"')

#-------------------------------------------------

if opts.verbose:
    print "reading in : %s"%inFITS
post = hp.read_map( inFITS, verbose=False )

#------------------------

if opts.source_coord=="C":
    if opts.verbose: 
        print "  interpretting as Equatorial coordinates (C)"

    if opts.target_coord=="C":
        if opts.verbose:
            print "  no rotation necessary (--target-coord=C)"
        
    elif opts.target_coord=="E":
        if opts.verbose:
            print "  rotating to Geographic coordinates (E) at time=%.3f"%gps
        post = triangulate.rotateMapC2E( post, gps )
        
    else:
        raise ValueError('--target-ccord must be either "C" or "E"')
        
elif opts.source_coord=="E":
    if opts.verbose:
        print "  interpretting as Geographic coordinates (E)"

    if opts.target_coord=="C":
        if opts.verbose:
            print "  rotating to Equatorial coordinates (C) at time=%.3f"%gps
        post = triangulate.rotateMapE2C( post, gps )

    elif opts.target_coord=="E":
        if opts.verbose:
            print "  no rotation necessary (--target-coord=E)"

    else:
        raise ValueError('--target-ccord must be either "C" or "E"')

else:
    raise ValueError('--source-coord must be either "C" or "E"')

#------------------------

if opts.verbose:
    print "writing to : %s"%outFITS
hp.write_map( outFITS, post )
