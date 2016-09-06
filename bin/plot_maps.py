#!/usr/bin/python
usage       = "compare_maps.py [--options] label1,fits1 label2,fits2 ..."
description = "plot skymaps on a figure for visualization. Basically a wrapper for lalinference.plot.healpix_heatmap"
author      = "R. Essick (reed.essick@ligo.org)"

#=================================================

import os
import json

import numpy as np
import healpy as hp

from plotting import mollweide as mw
plt = mw.plt
lalinf_plot = mw.lalinf_plot

from optparse import OptionParser

#================================================-

colors = ['b', 'r', 'g', 'm', 'c', 'k', 'y']

#=================================================

### parse arguments

parser = OptionParser(usage=usage, description=description)

parser.add_option("-v", "--verbose", default=False, action="store_true")

parser.add_option("", "--stack-posteriors", default=False, action="store_true")
parser.add_option("", "--stack-posteriors-background", default=None, type="string", help="a FITS file to plot in the background of the stacked plot")
parser.add_option("", "--stack-posteriors-linewidths", default=1, type="float", help="the linewidth for contours on stacked plot")
parser.add_option("", "--stack-posteriors-alpha", default=1.0, type="float", help="the alpha for contours on stacked plot")
parser.add_option("", "--stack-posteriors-levels", default=[], type='float', action='append', help='the confidence levels for contours in stacked plot')

parser.add_option("-l", "--logarithmic", default=False, action="store_true")

parser.add_option("-H", "--figheight", default=5, type="float")
parser.add_option("-W", "--figwidth", default=9, type="float")

parser.add_option("-o", "--output-dir", default=".", type="string")

parser.add_option("-p", "--projection", default="astro mollweide", type="string", help="either \"mollweide\", \"astro mollweide\"")
parser.add_option("", "--color-map", default="Reds", type="string", help="Default=\"Reds\"")

parser.add_option("-t", "--tag", default="", type="string")

parser.add_option("-T", "--transparent", default=False, action="store_true")

parser.add_option("", "--figtype", default=[], action="append", type="string")
parser.add_option("", "--dpi", default=500, type="int")

parser.add_option("", "--line-of-sight", default=[], action="append", type="string", help="eg: HL")
parser.add_option("", "--line-of-sight-color", default='k', type='string', help="the text and marker color for line-of-sight annotations")

parser.add_option("", "--zenith", default=[], action="append", type="string", help="eg: H")
parser.add_option("", "--zenith-color", default='k', type='string', help='the text and marker color for zenith annotations')

parser.add_option("", "--time-delay", default=[], action="append", type="string", help="eg: HL")
parser.add_option("", "--time-delay-Dec-RA", nargs=2, default=[], action="append", type="float", help="Should be specified in radians and this option requires two arguments (--time-delay-Dec-RA ${dec} ${ra}). If suppplied, we use this point to define time-delays (if told to plot them). If coord==C, this is interpreted as Dec,RA. If coord==E, this is interpreted as Theta,Phi")
parser.add_option("", "--time-delay-degrees", default=False, action="store_true", help="interpret --time-delay-Dec-RA as degrees")
parser.add_option("", "--time-delay-color", default='k', type='string', help='the line color for time-delay lines')
parser.add_option("", "--time-delay-alpha", default=1.0, type='float', help='the alpha saturation for time-delay lines')

parser.add_option("", "--marker-Dec-RA", nargs=2, default=[], action="append", type="float", help="Should be specified in adians and this option requires two arguments (--marker-Dec-RA ${dec} ${ra}). If suppplied, we label this point with a circles (if told to plot them). If coord==C, this is interpreted as Dec,RA. If coord==E, this is interpreted as Theta,Phi.")
parser.add_option("", "--marker-degrees", default=False, action="store_true", help="interpret --marker-Dec-RA as degrees")
parser.add_option("", "--marker-color", default='k', type='string', help='the edge-color for the markers')
parser.add_option("", "--marker-alpha", default=1.0, type='float', help='the alpha saturation for markers')
parser.add_option("", "--marker", default='o', type='string', help="the actual marker shape used")
parser.add_option("", "--marker-size", default=4, type='float', help='the size of the marker')
parser.add_option("", "--marker-edge-width", default=1, type='float', help='the edge width of the marker')

parser.add_option("", "--gps", default=None, type="float", help="must be specified if --line-of-sight or --zenith is used")
parser.add_option("", "--coord", default="C", type="string", help="coordinate system of the maps. Default is celestial (C), but we also know Earth-Fixed (E)")

parser.add_option("", "--continents", default=False, action="store_true", help="draw the continents on the map. Only used if --coord=E and --projection=mollweide")
parser.add_option("", "--continents-color", default='k', type='string', help='the color used to draw the continents')
parser.add_option("", "--continents-alpha", default=0.5, type='float', help='the alpha value for the contintents')

parser.add_option("", "--outline-labels", default=False, action="store_true", help="put a white outline around axis labels")

opts, args = parser.parse_args()

if not opts.figtype:
    opts.figtype.append( "png" )

if opts.tag:
    opts.tag = "_%s"%opts.tag

maps = {}
for arg in args:
    label, fits = arg.split(",")
    maps[label] = {"fits":fits}

labels = sorted(maps.keys())

if (opts.line_of_sight or opts.zenith or (opts.time_delay and opts.time_delay_Dec_RA)) and (opts.coord!="E") and (opts.gps==None):
    opts.gps = float(raw_input("gps = "))

opts.continents =  opts.continents and (opts.coord=="E") and (opts.projection=="mollweide")

if opts.stack_posteriors and (not opts.stack_posteriors_levels):
    opts.stack_posteriors_levels = [0.1, 0.5, 0.9]

#=============================================

### figure out positions for line-of-sight
line_of_sight = mw.gen_line_of_sight( opts.line_of_sight, coord=opts.coord, gps=opts.gps )

### figure out postions for zenith markers
zenith = mw.gen_zenith( opts.zenith, coord=opts.coord, gps=opts.gps )

### figure out points for time_delay
time_delay = mw.gen_time_delay( opts.time_delay_Dec_RA, opts.time_delay, coord=opts.coord, gps=opts.gps, degrees=opts.time_delay_degrees )

### figure out points for markers
marker_Dec_RA = mw.gen_marker_Dec_RA( opts.marker_Dec_RA, coord=opts.coord, gps=opts.gps, degrees=opts.marker_degrees )

#=============================================
### generate plots

figind = 0
if opts.stack_posteriors:
    stack_fig, stack_ax = mw.gen_fig_ax( figind, figwidth=opts.figwidth, figheight=opts.figheight, projection=opts.projection )
    figind += 1

    if opts.stack_posteriors_background:
        if opts.verbose:
            print "reading map from", fits
        post, header = hp.read_map( opts.stack_posteriors_background, h=True )
        if opts.verbose:
            print "plotting background for stackedPosteriors"
        mw.heatmap( post, stack_ax, color_map=opts.color_map )

    mw.annotate( stack_ax,
                 projection          = opts.projection,
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
                 continents          = opts.continents,
                 continents_color    = opts.continents_color,
                 continents_alpha    = opts.continents_alpha,
               )

for label in labels:
    fits = maps[label]['fits']
    if opts.verbose:
        print "reading map from", fits
    post, header = hp.read_map( fits, h=True )
    npix = len(post)
    if (dict(header)['ORDERING']=='NEST'): ### convert to RING ordering
        post = hp.nest2ring(nside, post)
    nside = hp.npix2nside(npix)
    if opts.verbose:
        print "\tnside=%d"%nside

    fig, ax = mw.gen_fig_ax( figind, figwidth=opts.figwidth, figheight=opts.figheight, projection=opts.projection )
    figind += 1

    mw.heatmap( post, ax, color_map=opts.color_map )
    mw.annotate( ax,
                 projection          = opts.projection,
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
                 continents          = opts.continents,
                 continents_color    = opts.continents_color,
                 continents_alpha    = opts.continents_alpha,
               )

    if opts.transparent:
        fig.patch.set_alpha(0.)
        ax.patch.set_alpha(0.)
        ax.set_alpha(0.)

    if opts.outline_labels:
        lalinf_plot.outline_text(ax)

    ### save the individual plot
    for figtype in opts.figtype:
        figname = "%s/%s%s.%s"%(opts.output_dir, label, opts.tag, figtype)
        if opts.verbose:
            print "\t", figname
        plt.savefig( figname, dpi=opts.dpi )
    plt.close( fig )

    ### plot contribution to stacked posteriors
    if opts.stack_posteriors:
        color = colors[(figind-1)%len(colors)]
        mw.contour( post, 
                    stack_ax, 
                    colors     = color, 
                    projection = opts.projection,
                    levels     = opts.stack_posterior_levels, 
                    alpha      = opts.stack_posteriors_alpha, 
                    linewidths = opts.stack_posteriors_linewidths )
        stack_fig.text(0.01, 0.99-0.05*(figind-1), label, color=color, ha='left', va='top')

### save the stacked posterior plot
if opts.stack_posteriors:
    plt.figure( 0 )
    plt.sca( stack_ax )

    if opts.transparent:
        stack_fig.patch.set_alpha(0.)
        stack_ax.patch.set_alpha(0.)
        stack_ax.set_alpha(0.)

    if opts.outline_labels:
        lalinf_plot.outline_text(stack_ax)

    for figtype in opts.figtype:
        figname = "%s/stackedPosterior%s.%s"%(opts.output_dir, opts.tag, figtype)
        if opts.verbose:
            print figname
        plt.savefig( figname, dpi=opts.dpi )
    plt.close( stack_fig )
