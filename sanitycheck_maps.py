description = """generate some sanity checks based on basic triangulation"""
usage = "sanitycheck_maps.py [--options] label,fits label,fits label,fits ..."
author = "Reed Essick (reed.essick@ligo.org)"

import triangulate
import healpy as hp
import numpy as np

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
plt.rcParams.update( {"text.usetex":True} )

from optparse import OptionParser

#=================================================

parser = OptionParser(usage=usage, description=description)

parser.add_option("-v", "--verbose", default=False, action="store_true")

parser.add_option("-L", "--line-of-sight", dest="los", defualt=[], action="append", type="string", help="ifo1,ifo2")

parser.add_option("-O", "--overhead", default=[], action="append", type="string", help="ifo")

parser.add_option("-o", "--output-dir", default=".", type="float")
parser.add_option("-t", "--tag", default="", type="float")

opts, args = parser.parse_args()

#=================================================

### read in maps
if opts.verbose:
    print "loading fits files:"
maps = {}
for arg in args:
    label, filename = arg.split(",")
    if opts.verbose:
        print "\t%s <- %s"%(label, filename)
    m = hp.read_map( filename )
    nside = hp.npix2nside( len(m) )
    if opts.verbose:
        print "\t\tnside = %d"%(nside)
    maps[label] = {"filename":filename, "map":m, "nside":nside}
labels = sorted(maps.keys())

#=================================================

### 1D line-of-sight
if opts.verbose:
    print "line-of-sight"
for opt in opts.los:
    ifo1, ifo2 = opt.split(",")
    if opts.verbose:
        print "\t %s -> %s"%(ifo1, ifo2)

    ### histograms of the distributions in coords defined by line-of-sight 
    ### make this like a corner plot, with projected histograms
    ### probably want to handle the formatting yourself, though
    ### also make a mollweide projection
"""
use this to look for whether the rings are "concentric" or "centered" where we expect them to be.
should show up as lines of constant lattitude
"""

#=================================================

### distance from overhead
if opts.verbose:
    print "overhead"
for ifo in opts.overhead:
    if opts.verbose:
        print "\t%s"%(ifo)

    ### histogram of the distributions in coords defined by overhead 
    ### again, make this a custom-made-corner plot.
    ### also make a mollweide projection
"""
use this to look for the "distance from the bulk" or other similar measures.
should allow us to measure how "unlikely" a certain sky location is
"""

