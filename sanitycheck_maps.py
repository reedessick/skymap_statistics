description = """generate some sanity checks based on basic triangulation"""
usage = "sanitycheck_maps.py [--options] label,fits label,fits label,fits ..."
author = "Reed Essick (reed.essick@ligo.org)"

import triangulate
import visualize
plt = visualize.plt

import healpy as hp
import numpy as np

from optparse import OptionParser

#=================================================

def compute_mi( theta, phi, Nbins, weights=None ):
    theta_bins = np.linspace(0, np.pi, Nbins+1)
    phi_bins = np.linspace(-np.pi, np.pi, Nbins+1)

    count = np.histogram2d( phi, theta, bins=(phi_bins, theta_bins), weights=weights )[0]

    return triangulate.mutualinformation( count, bins=None )
#    return triangulate.mutualinformation( count, bins=(theta_bins, phi_bins) )

#=================================================

colors = "b r g c m y k".split()

#=================================================

parser = OptionParser(usage=usage, description=description)

parser.add_option("-v", "--verbose", default=False, action="store_true")

parser.add_option("-L", "--line-of-sight", dest="los", default=[], action="append", type="string", help="ifo1,ifo2")
parser.add_option("-O", "--overhead", default=[], action="append", type="string", help="ifo")

parser.add_option("-p", "--plots", default=False, action="store_true", help="make plots in the new frame of reference")
parser.add_option("-m", "--mutualinformation", default=False, action="store_true", help="compute the mutual information in the new frame of reference")

parser.add_option("-c", "--coord", default="C", type="string", help="\"C\" or \"E\"")
parser.add_option("-T", "--t_geocent", default=None, type=float)

parser.add_option("-l", "--log", default=False, action="store_true", help="log scale histograms")
parser.add_option("-C", "--contour", default=False, action="store_true", help="plot contours instead of images")

parser.add_option("-g", "--grid", default=False, action="store_true")
parser.add_option("-o", "--output-dir", default=".", type="string")
parser.add_option("-t", "--tag", default="", type="string")

opts, args = parser.parse_args()

if opts.coord not in ["C", "E"]:
    raise ValueError("--coord=%s not understood"%(opts.coord))

if opts.coord=="C" and opts.t_geocent==None:
    opts.t_geocent = float( raw_input("t_geocent = ") )

if opts.tag:
    opts.tag = "_%s"%(opts.tag)

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
    npix = len(m)
    nside = hp.npix2nside( npix )
    if opts.verbose:
        print "\t\tnside = %d"%(nside)
    maps[label] = {"filename":filename, "map":m, "nside":nside, "npix":npix}

    if opts.mutualinformation:
        Nbins = max(100, int(npix**0.5/5))
        theta, phi = hp.pix2ang(nside, np.arange(npix))
        mi, entj = compute_mi( theta, phi, Nbins, weights=m )
        print "\tmutualinformationDistance(%s) : %.6f"%(label, mi/entj)
        maps[label].update( {'mutualinformation':mi/entj} )

labels = sorted(maps.keys())

#=================================================

### line-of-sight
if opts.verbose and len(opts.los):
    print "line-of-sight"
for opt in opts.los:
    ifo1, ifo2 = opt.split(",")
    if opts.verbose:
        print "\t %s -> %s"%(ifo1, ifo2)

    t, p = triangulate.line_of_sight( ifo1, ifo2, coord=opts.coord, tgeocent=opts.t_geocent )
    if opts.coord=="C":
        t = 0.5*np.pi - t ### convert from dec to theta    

    figax = None

    for cind, label in enumerate(labels):
        color = colors[cind%len(colors)]

        if opts.verbose:
            print "\t\t%s"%(label)
        m = maps[label]['map']
        nside = maps[label]['nside']
        npix = maps[label]['npix']
        
        theta, phi = hp.pix2ang(nside, np.arange(npix))
        ### rotate
        rtheta, rphi = triangulate.rotate2pole( theta, phi, t, p )

        Nbins = max(100, int(npix**0.5/5))

        if opts.mutualinformation:
            mi, entj = compute_mi( rtheta, rphi, Nbins, weights=m )
            print "\t\tmutualinformationDistance(%s) : %.6f"%(label, mi/entj)
        if opts.plots:
            figax = visualize.histogram2d( rtheta, rphi, Nbins=Nbins, weights=m, figax=figax, log=opts.log, contour=opts.contour, color=color )
            figax[0].text(0.99, 0.9-cind*0.05, label.replace('_','\_'), color=color, ha='right', va='top')

    if figax:
        fig, ax = figax
        ### decorate
	for a in ax:
		a.grid(opts.grid, which="both")
	plt.setp(ax[2].get_yticklabels(), visible=False)
	plt.setp(ax[1].get_xticklabels(), visible=False)

        figname = "%s/los-%s-%s%s.png"%(opts.output_dir, ifo1, ifo2, opts.tag)
        if opts.verbose:
            print "\t%s"%(figname)
        fig.savefig( figname )
        plt.close( fig )

"""
use this to look for whether the rings are "concentric" or "centered" where we expect them to be.
should show up as lines of constant lattitude
"""

#=================================================

### distance from overhead
if opts.verbose and len(opts.overhead):
    print "overhead"
for ifo in opts.overhead:
    if opts.verbose:
        print "\t%s"%(ifo)

    t, p = triangulate.overhead( ifo, coord=opts.coord, tgeocent=opts.t_geocent )
    if opts.coord=="C":
        t = 0.5*np.pi - t ### convert from dec to theta

    figax = None

    for cind, label in enumerate(labels):
        color = colors[cind%len(colors)]

        if opts.verbose:
            print "\t\t%s"%(label)
        m = maps[label]['map']
        nside = maps[label]['nside']
        npix = maps[label]['npix']

        theta, phi = hp.pix2ang(nside, np.arange(npix))
        ### rotate
        rtheta, rphi = triangulate.rotate2pole( theta, phi, t, p )

        Nbins = max(100, int(npix**0.5/5))

        if opts.mutualinformation:
            mi, entj = compute_mi( rtheta, rphi, Nbins, weights=m )
            print "\t\tmutualinformation(%s) : %.6f nats"%(label, mi/entj)
        if opts.plots:
            figax = visualize.histogram2d( rtheta, rphi, Nbins=Nbins, weights=m, figax=figax, log=opts.log, contour=opts.contour, color=color )
            figax[0].text(0.99, 0.9-cind*0.05, label.replace("_","\_"), color=color, ha='right', va='top')

    if figax:
        fig, ax = figax
        ### decorate
        for a in ax:
                a.grid(opts.grid, which="both")
        plt.setp(ax[2].get_yticklabels(), visible=False)
        plt.setp(ax[1].get_xticklabels(), visible=False)

        figname = "%s/overhead-%s%s.png"%(opts.output_dir, ifo, opts.tag)
        if opts.verbose:
            print "\t%s"%(figname)
        fig.savefig( figname )
        plt.close( fig )

"""
use this to look for the "distance from the bulk" or other similar measures.
should allow us to measure how "unlikely" a certain sky location is
"""

