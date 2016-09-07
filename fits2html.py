#!/usr/bin/python
usage       = "fits2html.py [--options] fits gps"
description = "generates plots and statistics summarizing a FITS file, and writes html documents to display the results"
author      = "reed.essick@ligo.org"

#-------------------------------------------------

import os
import json

import stats
import triangulate

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from plotting import mollweide as mw
from plotting import cartesian as ct

import detector_cache
import antenna

import numpy as np
import healpy as hp

from optparse import OptionParser

#-------------------------------------------------

parser = OptionParser(usage=usage, description=description)

#--- 

parser.add_option('-v', '--verbose', default=False, action='store_true')

# options about input FITS file

parser.add_option("", "--gps", default=None, type="float")
parser.add_option('-c', '--coord', default='C', type='string', help='either C or E. DEFAULT=C')

# options about network

parser.add_option('-i', '--ifo', default=[], type='string', action='append')

# general plotting options

parser.add_option('', '--color-map', default='OrRd', type='string', help='color map for heatmaps')
parser.add_option("-T", "--transparent", default=False, action="store_true")

parser.add_option("", "--time-delay-color", default='k', type='string', help='the line color for time-delay lines')
parser.add_option("", "--time-delay-alpha", default=1.0, type='float', help='the alpha saturation for time-delay lines')

# plotting options for mollweide

parser.add_option('', "--mollweide-levels", default=[], action='append', type='float', help='levels for mollweide countours')
parser.add_option('', '--mollweide-alpha', default=1.0, type='float', help='alpha for mollweide contours')
parser.add_option('', '--mollweide-linewidth', default=1.0, type='float', help='linewidth for mollweide contours')

parser.add_option("", "--line-of-sight-color", default='k', type='string', help="the text and marker color for line-of-sight annotations")

parser.add_option("", "--zenith-color", default='k', type='string', help='the text and marker color for zenith annotations')

parser.add_option("", "--marker-color", default='k', type='string', help='the edge-color for the markers')
parser.add_option("", "--marker-alpha", default=1.0, type='float', help='the alpha saturation for markers')
parser.add_option("", "--marker", default='o', type='string', help="the actual marker shape used")
parser.add_option("", "--marker-size", default=4, type='float', help='the size of the marker')
parser.add_option("", "--marker-edgewidth", default=1, type='float', help='the edge width of the marker')

# plotting options for dT marginals

parser.add_option("", "--dT-Nsamp", default=501, type='int')
parser.add_option("", "--dT-nside", default=None, type='int')
parser.add_option("", "--dT-xlim-dB", default=-20, type='float')
parser.add_option("", "--dT-no-yticks", default=False, action="store_true")

# options for computing statistics

parser.add_option("", "--base", default=2.0, type='float', help='base of logarithm used to compute entropy and information')
parser.add_option("", "--conf", default=[], action='append', type='float', help='the confidence levels used to evaluate size of map')

# general options

parser.add_option('-o', '--output-dir', default='.', type='string')
parser.add_option('-t', '--tag', default='', type='string')

parser.add_option("", "--figtype", default=[], action="append", type="string")
parser.add_option("", "--dpi", default=500, type="int")

#---

opts, args = parser.parse_args()

if opts.tag:
    opts.tag = "_"+opts.tag

if not opts.figtype:
    opts.figtype.append( "png" )

if not opts.mollweide_levels:
    opts.mollweide_levels = [0.10, 0.50, 0.90]

if not opts.conf:
    opts.conf = np.linspace(0, 1.0, 51)[1:]

for ifo in opts.ifo:
    if not detector_cache.detectors.has_key(ifo):
        raise ValueError("--ifo=%s not understood"%ifo)

if len(args)!=2:
    raise ValueError('please supply exactly 2 arguments\n%s'%usage)
fits = args[0]
gps  = float(args[1])

#-------------------------------------------------

opts.ifo = sorted(opts.ifo)
ifo_pairs = [ ifo1+ifo2 for ind, ifo1 in enumerate(opts.ifo) for ifo2 in opts.ifo[ind+1:] ]

#-------------------------------------------------

### create plots and summary info for map
outdir = os.path.join( opts.output_dir, os.path.basename(fits).split('.')[0] )
if not os.path.exists(outdir):
    os.makedirs(outdir)

if opts.verbose:
    print "%s -> %s"%(fits, outdir)

#-------------------------------------------------

data = {} ### store filenames and misc data
          ### probably want to replace this with some kind of class that knows how to write an html document and supporting data files

#-------------------------------------------------

### read in map and set up both coordinate systems
if opts.coord=="C":
    postC = hp.read_map( fits )
    postE = triangulate.rotateMapC2E( postC, gps )

elif opts.coord=="E":
    postE = hp.read_map( fits )
    postC = triangulate.rotateMapE2C( postE, gps )

else:
    raise ValueError('--coord=%s not understood'%opts.coord)

npix = len(postC)
nside = hp.npix2nside( npix )
theta, phi = hp.pix2ang( nside, np.arange(npix) )

#-------------------------------------------------

### generate plots

figind = 0

#-----------

if opts.verbose:
    print "generating mollweide projections"

for projection, coord, post in [('astro mollweide', 'C', postC), ('mollweide', 'E', postE)]:

    if opts.verbose:
        print "  cood=%s"%coord

    mapY = theta[mapIND]
    mapX = phi[mapIND]

    if coord == 'C':
        mapY = 0.5*np.pi - mapY ### convert from Theta->Dec
        X = triangulate.rotateRAC2E( phi, gps ) ### convert RA->phi
        Y = theta 
    else:
        X = phi
        Y = theta

    ### set up line-of-sight
    line_of_sight = mw.gen_line_of_sight( ifo_pairs, coord=coord, gps=gps )
                                          
    ### set up zenith
    zenith = mw.gen_zenith( opts.ifo, coord=coord, gps=gps )

    ### set up time-delay (for MAP)
    time_delay = mw.gen_time_delay( [(mapY, mapX)], ifo_pairs, coord=coord, gps=gps, degrees=False )

    ### set up marker (for MAP)
    marker_Dec_RA = mw.gen_marker_Dec_RA( [(mapY, mapX)], coord=coord, gps=gps, degrees=False )

    ### generate figure
    fig, ax = mw.gen_fig_ax( figind, projection='astro mollweide' )
    figind += 1

    mw.heatmap( post, ax, color_map=opts.color_map )
    mw.annotate( ax, 
                 continents       = coord=='E',
                 continents_color = opts.continents_color,
                 continents_alpha = opts.continents_alpha, )

    if opts.transparent:
        stack_fig.patch.set_alpha(0.)
        stack_ax.patch.set_alpha(0.)
        stack_ax.set_alpha(0.)

    ### save just the heatmap
    fignames = []
    for figtype in opts.figtype:
        figname = os.path.join(outdir, "heatmap%s%s.%s"%(coord, opts.tag, figtype) )
        if opts.verbose:
            print "    "+figname
        fig.savefig( figname, dpi=opts.dpi )
        fignames.append( figname )
    data['mw %s'%coord] = fignames
    
    ### annotate with fancy crap
    mw.annoate( ax, 
                projection          = projection,
                line_of_sight       = line_of_sight,
                line_of_sight_color = opts.line_of_sight_color,
                zenith              = zenith,
                zenith_color        = opts.zenith_color,
                time_delay          = time_delay,
                time_delay_color    = opts.time_delay_color,
                time_delay_alpha    = opts.time_delay_alpha,
                marker_Dec_RA       = marker_Dec_RA,
                marker              = opts.marker,
                marker_color        = opts.marker_color,
                marker_size         = opts.marker_size,
                marker_edgewidth    = opts.marker_edgewidth,
                marker_alpha        = opts.marker_alpha,
              )

    ### save heatmap + fancy crap
    fignames = []
    for figtype in opts.figtype:
        figname = os.path.join(outdir, "heatmap%s-annotated%s.%s"%(coord, opts.tag, figtype) )
        if opts.verbose:
            print "    "+figname
        fig.savefig( figname, dpi=opts.dpi )
        fignames.append( figname )
    data['mw %s ann'%coord] = fignames

    ### add antenna patterns as contours
    for ifo in opts.ifo:
        if not detector_cache.detectors.has_key(ifo):
            raise ValueError("ifo=%s not understood"%(ifo))

        Fp, Fx = detector_cache.detectors[ifo].antenna_patterns( Y, X, 0.0 )
        ant = Fp**2 + Fx**2
        ant /= np.sum(ant)

        mw.contour( ant,
                    ax,
                    colors     = color,
                    levels     = opts.mollweide_levels,
                    alpha      = opts.mollweide_alpha,
                    linewidths = opts.mollweide_linewidths )
        fig.text(0.01, 0.99-0.05*(figind-1), ifo, color=color, ha='left', va='top')

    ### save heatmap + fancy crap + antenna pattern contours
    fignames = []
    for figtype in opts.figtype:
        figname = os.path.join(outdir, "heatmap%s-antennas%s.%s"%(coord, opts.tag, figtype) )
        if opts.verbose:
            print "    "+figname
        fig.savefig( figname, dpi=opts.dpi )
        fignames.append( figname )
    data['mw %s ant'%coord] = fignames

    ### done with this figure
    plt.close( fig )

#-----------

if opts.verbose:
    print "generating dT marginals"

for ifo1, ifo2 in pairs:
    ifos = "".join([ifo1, ifo2])
    if opts.verbose:
        print "  %s - %s"%(ifo1, ifo2)

    sampDt = ct.get_sampDt( ifos, Nsamp=opts.dT_Nsamp )
    maxDt = sampDt[-1]

    fig, ax = ct.genDT_fig_ax( figind )
    figind += 1

    ax.set_xlim(xmin=maxDt*1e3, xmax=-maxDt*1e3) ### we work in ms here...

    if opts.dT_nside:
        kde = ct.post2marg( stats.resample(postE, opts.dT_nside), ifos, sampDt, coord='E', gps=opts.gps )
    else:
        kde = ct.post2marg( postE, ifos, sampDt, coord='E', gps=opts.gps )

    ### plot
    ct.plot( ax, sampDt, kde, xlim_dB=opts.dT_xlim_dB )

    ### decorate
    ax.set_xlabel(r'$\Delta t_{%s}\ [\mathrm{ms}]$'%(ifos))
    ax.set_ylabel(r'$p(\Delta t_{%s}|\mathrm{data})$'%(ifos))

    if opts.dT_no_yticks:
        ax.set_yticklabels([])

    if opts.transparent:
        stack_fig.patch.set_alpha(0.)
        stack_ax.patch.set_alpha(0.)
        stack_ax.set_alpha(0.)

    ### save just dT marginals
    fignames = []
    for figtype in opts.figtype:
        figname = os.path.join(outdir, "dT_%s%s.%s"%(ifos, opts.tag, figtype) )
        if opts.verbose:
            print figname
        fig.savefig(figname, dpi=opts.dpi)
        fignames.append( figname )
    data['dT %s'%ifos] = fignames

    ### annotate the plot
    ct.annotate( ax,
                 [ hp.pix2ang( nside, np.argmax(postE) ) ],
                 ifos,
                 maxDt,
                 coord   = 'E',
                 gps     = opts.gps,
                 color   = opts.time_delay_color,
                 alpha   = opts.time_delay_alpha,
                 degrees = opts.False,
               )

    ### save annotated dT marginals
    fignames = []
    for figtype in opts.figtype:
        figname = os.path.join(outdir, "dT_%s-annotated%s.%s"%(ifos, opts.tag, figtype) )
        if opts.verbose:
            print figname
        fig.savefig(figname, dpi=opts.dpi)
        fignames.append( figname )
    data['dT %s ann'%ifos] = fignames

    plt.close(fig)

#-----------

if opts.verbose:
    print "generating cartesian projections in line-of-sight frame"

for ifo1, ifo2 in pairs:
    if opts.verbose:
        print "  %s - %s"%(ifo1, ifo2)

    t, p = triangulate.line_of_sight( ifo1, ifo2, coord=opts.coord, tgeocent=gps )
    if opts.coord=="C":
        t = 0.5*np.pi - t

    fig, ax, rproj, tproj = ct.genHist_fig_ax( figind, figwidth=opts.figwidth, figheight=opts.figheight )
    figind += 1

    ### rotate
    rtheta, rphi = triangulate.rotate2pole( theta, phi, t, p )

    Nbins = max(100, int(npix**0.5/5))

    ### compute mutual info
    mi, entj = triangulate.compute_mi( rtheta, rphi, Nbins, weights=m )
    data['zen %s%s mi'%(ifo1,ifo2)] = mi
    data['zen %s%s Hj'%(ifo1,ifo2)] = entj

    ### plot
    ct.histogram2d( rtheta, rphi, ax, rproj, tproj, Nbins=Nbins, weights=m, log=opts.log, contour=opts.contour, color=color, cmap=opts.color_map )
    fig.text(0.99, 0.9-cind*0.05, label.replace('_','\_'), color=color, ha='right', va='top')

    ### decorate
    plt.setp(rproj.get_yticklabels(), visible=False)
    plt.setp(tproj.get_xticklabels(), visible=False)

    ### save
    fignames = []
    for figtype in opts.figtype:
        figname = "%s/los-%s-%s%s.%s"%(opts.output_dir, ifo1, ifo2, opts.tag, figtype)
        if opts.verbose:
            print "\t%s"%(figname)
        fig.savefig( figname )
        fignames.append( figname )
    data['los %s%s'%(ifo1,ifo2)] = fignames

    plt.close( fig )

#-----------

if opts.verbose:
    print "generating cartesian projections in overhead frame"

for ifo in opts.ifo:
    if opts.verbose:
        print "  %s"%ifo

    t, p = triangulate.overhead( ifo, coord=opts.coord, tgeocent=gps )
    if opts.coord=="C":
        t = 0.5*np.pi - t

    fig, ax, rproj, tproj = ct.genHist_fig_ax( figind, figwidth=opts.figwidth, figheight=opts.figheight )
    figind += 1

    ### rotate
    rtheta, rphi = triangulate.rotate2pole( theta, phi, t, p )

    Nbins = max(100, int(npix**0.5/5))

    ### compute mutual info
    mi, entj = triangulate.compute_mi( rtheta, rphi, Nbins, weights=m )
    data['zen %s mi'%ifo] = mi
    data['zen %s Hj'%ifo] = entj

    ### plot
    ct.histogram2d( rtheta, rphi, ax, rproj, tproj, Nbins=Nbins, weights=m, log=opts.log, contour=opts.contour, color=color, cmap=opts.color_map )
    fig.text(0.99, 0.9-cind*0.05, label.replace('_','\_'), color=color, ha='right', va='top')

    ### decorate
    plt.setp(rproj.get_yticklabels(), visible=False)
    plt.setp(tproj.get_xticklabels(), visible=False)

    ### save
    fignames = []
    for figtype in opts.figtype:
        figname = "%s/los-%s-%s%s.%s"%(opts.output_dir, ifo1, ifo2, opts.tag, figtype)
        if opts.verbose:
            print "\t%s"%(figname)
        fig.savefig( figname )
        fignames.append( figname )
    data['zen %s'%ifo] = fignames

    plt.close( fig )

#-------------------------------------------------

### compute statistics

if opts.verbose:
    print "computing statistics"

data['nside'] = nside

### information measures
if opts.verbose:
    print "  information measures"
data['H'] = stats.entropy( postC, base=opts.base )
data['I'] = stats.information( postC, base=opts.base )

### size of confidence regions
if opts.verbose:
    print "  confidence regions"
pixarea = hp.nside2pixarea( nside )
cr = stats.credible_region(postC, opts.conf)
cr_size = [len(_)*pixarea for _ in cr]
data['CR size'] = zip( opts.conf, cr_size )

raise NotImplementedError('write plot for confidence region size')

### antenna patterns
if opts.verbose:
    print "  antenna patterns"

mapind = np.argmax(postE)
for ifo in opts.ifo:
    if opts.verbose:
        print "    %s"%ifo

    Fp, Fx = detector_cache.detectors[ifo].antenna_patterns( theta, phi, 0.0 )
    
    ### MAP values
    raise NotImplementedError('compute values of Fp, Fx, Fp^2+Fx^2 at MAP')

    ### averaged over posterior
    raise NotImplementedError('compute values of Fp, Fx, Fp^2+Fx^2 averaged over the entire sky with posterior')

#-------------------------------------------------

### if we have posterior_samples.dat, generate postviz (only relevant for LIB and LALInf)

raise NotImplementedError('check for posterior_samples.dat and generate postviz if present')

#-------------------------------------------------

### write HTML document

raise NotImplementedError('write html document and supporting json data')

'''
want to make "data" a class with some method that does this for us cleanly.
'''
