#!/usr/bin/python
usage       = "fits2html.py [--options] fits gps"
description = "writes html documents summarizing fits files"
author      = "R. Essick (reed.essick@ligo.org)"

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

parser.add_option('-v', '--verbose', default=False, action='store_true')

parser.add_option('-c', '--coord', default='C', type='string', help='either C or E. DEFAULT=C')

parser.add_option('-i', '--ifo', default=[], type='string', action='append')

parser.add_option('', '--color-map', default='OrRd', type='string', help='color map for heatmaps')

parser.add_option("", "--line-of-sight-color", default='k', type='string', help="the text and marker color for line-of-sight annotations")

parser.add_option("", "--zenith-color", default='k', type='string', help='the text and marker color for zenith annotations')

parser.add_option("", "--time-delay-color", default='k', type='string', help='the line color for time-delay lines')
parser.add_option("", "--time-delay-alpha", default=1.0, type='float', help='the alpha saturation for time-delay lines')

parser.add_option("", "--marker-color", default='k', type='string', help='the edge-color for the markers')
parser.add_option("", "--marker-alpha", default=1.0, type='float', help='the alpha saturation for markers')
parser.add_option("", "--marker", default='o', type='string', help="the actual marker shape used")
parser.add_option("", "--marker-size", default=4, type='float', help='the size of the marker')
parser.add_option("", "--marker-edge-width", default=1, type='float', help='the edge width of the marker')

parser.add_option('-o', '--output-dir', default='.', type='string')
parser.add_option('-t', '--tag', default='', type='string')

parser.add_option("", "--figtype", default=[], action="append", type="string")
parser.add_option("", "--dpi", default=500, type="int")

opts, args = parser.parse_args()

if opts.tag:
    opts.tag = "_"+opts.tag

if not opts.figtype:
    opts.figtype.append( "png" )

if len(args)!=2:
    raise ValueError('please supply exactly 2 arguments\n%s'%usage)
fits = args[0]
gps  = float(args[1])

#-------------------------------------------------

ifo_pairs = [ ifo1+ifo2 for ind, ifo1 in enumerate(opts.ifo) for ifo2 in opts.ifo[ind+1:] ]

#-------------------------------------------------

raise NotImplementedError('WRITE ME!')








### create plots and summary info for map
outdir = os.path.join( opts.output_dir, os.path.basename(fits).split('.')[0] )
if not os.path.exists(outdir):
    os.makedirs(outdir)

if opts.verbose:
    print "%s -> %s"%(fits, outdir)

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

### generate plots
figing = 0
if opts.verbose:
    print "generating mollweide projections"

### Equatorial projections
for projection, coord, post in [('astro mollweide', 'C', postC), ('mollweide', 'E', postE)]:

    npix = len(post)
    nside = hp.npix2nside( npix )
    theta, phi = hp.pix2ang( nside, np.arange(npix) )

    mapY = theta[mapIND]
    mapX = phi[mapIND]

    if coord == 'C':
        mapY = 0.5*np.pi - mapY ### convert from Theta->Dec
        phi = triangulate.rotateRAC2E( phi, gps ) ### convert RA->phi

    ### set up line-of-sight
    line_of_sight = pm.gen_line_of_sight( ifo_pairs, coord=coord, gps=gps )
                                          
    ### set up zenith
    zenith = pm.gen_zenith( opts.ifo, coord=coord, gps=gps )

    ### set up time-delay (for MAP)
    time_delay = gen_time_delay( [(mapY, mapX)], ifo_pairs, coord=coord, gps=gps, degrees=False )

    ### set up marker (for MAP)
    marker_Dec_RA = gen_marker_Dec_RA( [(mapY, mapX)], coord=coord, gps=gps, degrees=False )

    ### generate figure
    fig, ax = pm.gen_figax( figind, projection='astro mollweide' )
    figind += 1

    pm.heatmap( post, ax, color_map=opts.color_map )
    pm.annotate( ax, 
                 continents       = coord=='E',
                 continents_color = opts.continents_color,
                 continents_alpha = opts.continents_alpha, )

    ### save just the heatmap
    for figtype in opts.figtype:
        figname = os.path.join(outdir, "heatmapC%s.%s"%(opts.tag, figtype)
        if opts.verbose:
            print "    "+figname
        fig.savefig( figname, dpi=opts.dpi )
    
    ### annotate with fancy crap
    pm.annoate( ax, 
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
    for figtype in opts.figtype:
        figname = os.path.join(outdir, "heatmapC-annotated%s.%s"%(opts.tag, figtype)
        if opts.verbose:
            print "    "+figname
        fig.savefig( figname, dpi=opts.dpi )

    ### add antenna patterns as contours

    for ifo in opts.ifo:
        if not aa.detector_cache.detectors.has_key(ifo):
            raise ValueError("ifo=%s not understood"%(ifo))

        Fp, Fx = aa.detector_cache.detectors[ifo].antenna_patterns( theta, phi, 0.0 )
        ant = Fp**2 + Fx**2
        ant /= np.sum(ant)

        contour( ant,
                 ax,
                 colors     = color,
                 projection = projection,
                 levels     = opts.stack_posterior_levels,
                 alpha      = opts.stack_posteriors_alpha,
                 linewidths = opts.stack_posteriors_linewidths )
        fig.text(0.01, 0.99-0.05*(figind-1), ifo, color=color, ha='left', va='top')

    ### save heatmap + fancy crap + antenna pattern contours
    for figtype in opts.figtype:
        figname = os.path.join(outdir, "heatmapC-antennas%s.%s"%(opts.tag, figtype)
        if opts.verbose:
            print "    "+figname
        fig.savefig( figname, dpi=opts.dpi )

    ### done with this figure
    plt.close( fig )
















    '''
dT marginals

cartesian sanityCheck in line-of-sight frame
    '''

    ### compute statistics
    '''
nside
entropy
size of confidence regions (plot conf vs deg2 in addition to a table?)
mutual information distance in line-of-sight frame
estimates of how much skymap is accessible from different locations with different openning angles (php form to allow interactivity?)
values of Fp, Fx, Fp^2+Fx^2 at MAP, averaged over the entire sky with posterior
    '''

    ### if we have posterior_samples.dat, generate postviz
'''
only relevant for LIB and LALInf
'''
