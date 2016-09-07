#!/usr/bin/python
usage       = "sanitycheck_maps.py [--options] label,fits label,fits label,fits ..."
description = "generate some sanity checks based on basic triangulation"
author      = "reed.essick@ligo.org"

#=================================================

import triangulate

import healpy as hp
import numpy as np

from plotting import cartesian as ct
plt = ct.plt

from optparse import OptionParser

#=================================================

colors = "b r g c m y k".split()

#=================================================

parser = OptionParser(usage=usage, description=description)

parser.add_option("-v", "--verbose", default=False, action="store_true")

parser.add_option("-H", "--figheight", default=5, type="float")
parser.add_option("-W", "--figwidth", default=9, type="float")

parser.add_option("-L", "--line-of-sight", dest="los", default=[], action="append", type="string", help="eg: HL")
parser.add_option("-O", "--zenith", default=[], action="append", type="string", help="eg: H")

parser.add_option("-p", "--plots", default=False, action="store_true", help="make plots in the new frame of reference")
parser.add_option("-m", "--mutualinformation", default=False, action="store_true", help="compute the mutual information in the new frame of reference")

parser.add_option("-c", "--coord", default="C", type="string", help="\"C\" or \"E\"")
parser.add_option("-T", "--gps", default=None, type=float)

parser.add_option("-l", "--log", default=False, action="store_true", help="log scale histograms")
parser.add_option("-C", "--contour", default=False, action="store_true", help="plot contours instead of images")

parser.add_option("-o", "--output-dir", default=".", type="string")
parser.add_option("-t", "--tag", default="", type="string")

parser.add_option("", "--figtype", default=[], action="append", type="string")
parser.add_option("", "--dpi", default=500, type="int")

parser.add_option("", "--color-map", default="Reds", type="string", help="Default=\"Reds\"")

opts, args = parser.parse_args()

if opts.coord not in ["C", "E"]:
    raise ValueError("--coord=%s not understood"%(opts.coord))

if opts.coord=="C" and opts.gps==None:
    opts.gps = float( raw_input("gps = ") )

if opts.tag:
    opts.tag = "_%s"%(opts.tag)

if not opts.figtype:
    opts.figtype.append( "png" )

#=================================================

### read in maps
if opts.verbose:
    print "loading fits files:"
maps = {}
for arg in args:
    label, filename = arg.split(",")
    if opts.verbose:
        print "\t%s <- %s"%(label, filename)
    m = hp.read_map( filename, verbose=False )
    npix = len(m)
    nside = hp.npix2nside( npix )
    if opts.verbose:
        print "\t\tnside = %d"%(nside)
    maps[label] = {"filename":filename, "map":m, "nside":nside, "npix":npix}

    if opts.mutualinformation:
        Nbins = max(100, int(npix**0.5/5))
        theta, phi = hp.pix2ang(nside, np.arange(npix))
        mi, entj = triangulate.compute_mi( theta, phi, Nbins, weights=m )
        print "\tmutualinformationDistance(%s) : %.6f"%(label, mi/entj)
        maps[label].update( {'mutualinformation':mi/entj} )

labels = sorted(maps.keys())

#=================================================

figind = 0

### line-of-sight
if opts.verbose and len(opts.los):
    print "line-of-sight"
for opt in opts.los:
    ifo1, ifo2 = opt
    if opts.verbose:
        print "\t %s -> %s"%(ifo1, ifo2)

    t, p = triangulate.line_of_sight( ifo1, ifo2, coord=opts.coord, tgeocent=opts.gps )
    if opts.coord=="C":
        t = 0.5*np.pi - t ### convert from dec to theta    

    fig, ax, rproj, tproj = ct.genHist_fig_ax( figind, figwidth=opts.figwidth, figheight=opts.figheight )
    figind += 1

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
            mi, entj = triangulate.compute_mi( rtheta, rphi, Nbins, weights=m )
            print "\t\tmutualinformationDistance(%s) : %.6f"%(label, mi/entj)

        if opts.plots:
            ct.histogram2d( rtheta, rphi, ax, rproj, tproj, Nbins=Nbins, weights=m, log=opts.log, contour=opts.contour, color=color, cmap=opts.color_map )
            fig.text(0.99, 0.9-cind*0.05, label.replace('_','\_'), color=color, ha='right', va='top')

    ### decorate
    plt.setp(rproj.get_yticklabels(), visible=False)
    plt.setp(tproj.get_xticklabels(), visible=False)

    for figtype in opts.figtype:
        figname = "%s/los-%s-%s%s.%s"%(opts.output_dir, ifo1, ifo2, opts.tag, figtype)
        if opts.verbose:
            print "\t%s"%(figname)
        fig.savefig( figname, dpi=opts.dpi )
    plt.close( fig )

#=================================================

### distance from zenith
if opts.verbose and len(opts.zenith):
    print "zenith"
for ifo in opts.zenith:
    if opts.verbose:
        print "\t%s"%(ifo)

    t, p = triangulate.overhead( ifo, coord=opts.coord, tgeocent=opts.gps )
    if opts.coord=="C":
        t = 0.5*np.pi - t ### convert from dec to theta

    fig, ax, rproj, tproj = ct.genHist_fig_ax( figind, figwidth=opts.figwidth, figheight=opts.figheight )
    figind += 1

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
            mi, entj = triangulate.compute_mi( rtheta, rphi, Nbins, weights=m )
            print "\t\tmutualinformation(%s) : %.6f nats"%(label, mi/entj)
        if opts.plots:
            ct.histogram2d( rtheta, rphi, ax, rproj, tproj, Nbins=Nbins, weights=m, log=opts.log, contour=opts.contour, color=color, cmap=opts.colormap )
            fig.text(0.99, 0.9-cind*0.05, label.replace("_","\_"), color=color, ha='right', va='top')

    ### decorate
    plt.setp(rproj.get_yticklabels(), visible=False)
    plt.setp(tproj.get_xticklabels(), visible=False)

    for figtype in opts.figtype:
        figname = "%s/zenith-%s%s.%s"%(opts.output_dir, ifo, opts.tag, figtype)
        if opts.verbose:
            print "\t%s"%(figname)
        fig.savefig( figname, dpi=opts.dpi )
    plt.close( fig )
