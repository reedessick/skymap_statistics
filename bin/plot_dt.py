#!/usr/bin/python
usage       = "plot_dt.py [--options] label,FITS label,FITS ..."
description = "plot the marigals for time of flight between the two specified IFOs"
author      = "R. Essick (reed.essick@ligo.org)"

#-------------------------------------------------

import os

import stats
import triangulate 

import healpy as hp
import numpy as np

from plotting import cartesian as ct
plt = ct.plt
from plotting import colors

from optparse import OptionParser

#-------------------------------------------------

### parse arguments

parser = OptionParser(usage=usage, description=description)

parser.add_option("-v", "--verbose", default=False, action='store_true')

parser.add_option("", "--stack-posteriors", default=False, action="store_true")

parser.add_option("-H", "--figheight", default=5, type="float")
parser.add_option("-W", "--figwidth", default=9, type="float")

parser.add_option("-o", "--output-dir", default=".", type="string")
parser.add_option("-t", "--tag", default="", type="string")

parser.add_option("-T", "--transparent", default=False, action="store_true")

parser.add_option("", "--time-delay", default=[], action="append", type="string", help="eg: HL")

parser.add_option("", "--time-delay-Dec-RA", nargs=2, default=[], action="append", type="float", help="Should be specified in radians and this option requires two arguments (--time-delay-Dec-RA ${dec} ${ra}). If suppplied, we use this point to define time-delays (if told to plot them). If coord==C, this is interpreted as Dec,RA. If coord==E, this is interpreted as Theta,Phi")
parser.add_option("", "--time-delay-degrees", default=False, action="store_true", help="interpret --time-delay-Dec-RA as degrees")
parser.add_option("", "--time-delay-color", default='k', type='string', help='the line color for time-delay lines')
parser.add_option("", "--time-delay-alpha", default=1.0, type='float', help='the alpha saturation for time-delay lines')

parser.add_option("", "--gps", default=None, type="float", help="must be specified if --line-of-sight or --zenith is used")
parser.add_option("", "--coord", default="C", type="string", help="coordinate system of the maps. Default is celestial (C), but we also know Earth-Fixed (E)")

parser.add_option("", "--figtype", default=[], action="append", type="string")
parser.add_option("", "--dpi", default=500, type="int")

parser.add_option("", "--Nsamp", default=501, type='int')
parser.add_option("", "--nside", default=None, type='int')
parser.add_option("", "--xlim-dB", default=-20, type='float')
parser.add_option("", "--no-yticks", default=False, action="store_true")

opts, args = parser.parse_args()

if not opts.time_delay:
    raise ValueError('please supply at least one IFO pair via --time-delay\n%s'%usage)

if not opts.figtype:
    opts.figtype.append( "png" )

if opts.tag:
    opts.tag = "_%s"%opts.tag

if not os.path.exists(opts.output_dir):
    os.makedirs(opts.output_dir)

if (opts.coord=="C") and (opts.gps==None):
    opts.gps = float(raw_input("gps = "))

#---------------------------------------------

if opts.verbose:
    print "reading maps"
maps = {}
for arg in args:
    label, fits = arg.split(",")
    if opts.verbose:
        print "    %s -> %s"%(fits, label)
    post = hp.read_map(fits, verbose=False)
    if opts.nside!=None:
        post = stats.resample( post, opts.nside )
    maps[label] = {"fits"  : fits, 
                   "post"  : post, 
                  }
labels = sorted(maps.keys())

#---------------------------------------------

figind = 0
for ifos in opts.time_delay:
    if opts.verbose:
        print "time delay between %s and %s"%(ifos[0], ifos[1])

    ### get max dt allowed
    sampDt = ct.gen_sampDt( ifos, Nsamp=opts.Nsamp )
    maxDt = sampDt[-1]

    ### compute and plot marginals
    if opts.stack_posteriors:
        stack_fig, stack_ax = ct.genDT_fig_ax( figind, figwidth=opts.figwidth, figheight=opts.figheight )
        figind += 1

        genColor = colors.getColor()

        stack_ax.set_xlim(xmin=maxDt*1e3, xmax=-maxDt*1e3) ### we work in ms here...

    for cind, label in enumerate(labels):
        if opts.verbose:
            print "    "+label

        post  = maps[label]['post']

        kde = ct.post2marg( post, ifos, sampDt, coord=opts.coord, gps=opts.gps ) 

        ### plot
        fig, ax = ct.genDT_fig_ax( figind, figwidth=opts.figwidth, figheight=opts.figheight )
        figind += 1

        ax.set_xlim(xmin=maxDt*1e3, xmax=-maxDt*1e3) ### we work in ms here...

        ct.plot( ax, sampDt, kde, label=label, color='b', xlim_dB=opts.xlim_dB )

        ### annotate the plot
        ct.annotate( ax, 
                     opts.time_delay_Dec_RA, 
                     ifos,
                     maxDt,  
                     coord   = opts.coord, 
                     gps     = opts.gps, 
                     color   = opts.time_delay_color, 
                     alpha   = opts.time_delay_alpha, 
                     degrees = opts.time_delay_degrees,
                   )

        ### decorate
        ax.set_xlabel(r'$\Delta t_{%s}\ [\mathrm{ms}]$'%(ifos))
        ax.set_ylabel(r'$p(\Delta t_{%s}|\mathrm{data})$'%(ifos))

        if opts.no_yticks:
            ax.set_yticklabels([])

        if opts.transparent:
            fig.patch.set_alpha(0.)
            ax.patch.set_alpha(0.)
            ax.set_alpha(0.)

        ### save figure
        for figtype in opts.figtype:
            figname = "%s/dT-%s_%s%s.%s"%(opts.output_dir, ifos, label, opts.tag, figtype)
            if opts.verbose:
                print "      "+figname
            fig.savefig(figname, dpi=opts.dpi)
        plt.close(fig)

        if opts.stack_posteriors:
            color = genColor.next()
            ct.plot( stack_ax, sampDt, kde, label=label, color=color, xlim_dB=opts.xlim_dB )
            stack_fig.text(0.10+0.02, 0.93-0.05*cind, label, color=color, ha='left', va='top')

    if opts.stack_posteriors:
        ### annotate the plot
        ct.annotate( stack_ax, 
                     opts.time_delay_Dec_RA, 
                     ifos,
                     maxDt,  
                     coord   = opts.coord, 
                     gps     = opts.gps, 
                     color   = opts.time_delay_color, 
                     alpha   = opts.time_delay_alpha, 
                     degrees = opts.time_delay_degrees,
                   )

        ### decorate
        stack_ax.set_xlabel(r'$\Delta t_{%s}\ [\mathrm{ms}]$'%(ifos))
        stack_ax.set_ylabel(r'$p(\Delta t_{%s}|\mathrm{data})$'%(ifos))

        if opts.no_yticks:
            stack_ax.set_yticklabels([])

        if opts.transparent:
            stack_fig.patch.set_alpha(0.)
            stack_ax.patch.set_alpha(0.)
            stack_ax.set_alpha(0.)

        ### save figure
        for figtype in opts.figtype:
            figname = "%s/dT-%s_stacked%s.%s"%(opts.output_dir, ifos, opts.tag, figtype)
            if opts.verbose:
                print figname
            stack_fig.savefig(figname, dpi=opts.dpi)
        plt.close(stack_fig)
