#!/usr/bin/python
usage = "plot_dTMarginals.py [--options] label,FITS label,FITS ..."
description = "plot the marigals for time of flight between the two specified IFOs"
author = "reed.essick@ligo.org"

#-------------------------------------------------

import stats
import triangulate 

import healpy as hp
import numpy as np

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
plt.rcParams.update({
  'font.family' : 'serif',
   })

from optparse import OptionParser

#-------------------------------------------------

twopi = 2*np.pi
pi2 = 0.5*np.pi

colors = ['b', 'r', 'g', 'm', 'c', 'k', 'y']

axpos = [0.10, 0.10, 0.85, 0.85]

#-------------------------------------------------

parser = OptionParser(usage=usage, description=description)

parser.add_option("-v", "--verbose", default=False, action='store_true')

parser.add_option("-H", "--figheight", default=5, type="float")
parser.add_option("-W", "--figwidth", default=9, type="float")

parser.add_option("-o", "--output-dir", default=".", type="string")
parser.add_option("-t", "--tag", default="", type="string")

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

if (opts.coord=="C") and (opts.gps==None):
    opts.gps = float(raw_input("gps = "))

#-------------------------------------------------

if opts.verbose:
    print "reading maps"
maps = {}
for arg in args:
    label, fits = arg.split(",")
    post = hp.read_map(fits, verbose=False)
    if opts.nside!=None:
        post = stats.resample( post, opts.nside )
        nside = opts.nside
    else:
        nside = hp.npix2nside(len(post))
    dAng = hp.nside2pixarea(nside, degrees=False)**0.5

    maps[label] = {"fits"  : fits, 
                   "post"  : post, 
                   "nside" : nside,
                   "dAng"  : dAng,
                  }

labels = sorted(maps.keys())

#-------------------------------------------------

for ifos in opts.time_delay:
    if opts.verbose:
        print "time delay between %s and %s"%(ifos[0], ifos[1])

    fig = plt.figure(figsize=(opts.figwidth, opts.figheight))
    ax = fig.add_axes(axpos)

    ### get max dt allowed
    t, p = triangulate.line_of_sight( ifos[0], ifos[1], coord='E' )
    maxDt = np.abs( triangulate.time_delay( t, p, ifos[0], ifos[1], coord='E' ) )
    sampDt = np.linspace(-maxDt, maxDt, opts.Nsamp)
    spacer = np.ones_like(sampDt)
    kde = np.zeros_like(sampDt)

    xmin = maxDt
    xmax = -maxDt

    cind = 0
    for label in labels:
        if opts.verbose:
            print "    ", label

        post  = maps[label]['post']
        nside = maps[label]['nside']
        dAng  = maps[label]['dAng']

        y, x = hp.pix2ang(nside, np.arange(len(post)))
        if opts.coord=="C":
            y = pi2 - y

        ### compute time-of-flight
        dt = triangulate.time_delay( y, x, ifos[1], ifos[0], coord=opts.coord, tgeocent=opts.gps, degrees=False )

        ### build up KDE
        sigma = maxDt * ( 1 - dt/maxDt )**0.5 * dAng

        big_post   = np.outer( spacer, post )
        big_dt     = np.outer( spacer, dt )
        big_sigma  = np.outer( spacer, sigma )
        big_sampDt = np.outer( sampDt, np.ones_like(dt) )

        kde = np.sum( big_post * np.exp( -( (big_dt - big_sampDt)/big_sigma )**2 ) / (twopi**0.5 * big_sigma), axis=-1 )
        kde /= np.sum(kde)

        subset = sampDt[kde>10**(opts.xlim_dB/10)*np.max(kde)]
        xmax = max(xmax, np.max(subset))
        xmin = min(xmin, np.min(subset))

        ### plot
        color = colors[cind%len(colors)]
        ax.plot( sampDt*1e3, kde*1e3, label=label, color=color )

        fig.text( axpos[0]+0.01*axpos[2], axpos[1]+0.99*axpos[-1]-cind*0.05, label, color=color, ha='left', va='top' )

        cind += 1

    ### add markers
    for dec, ra in opts.time_delay_Dec_RA:
        if opts.time_delay_degrees:
            dec *= triangulate.deg2rad
            ra *= triangulate.deg2rad
        dt = triangulate.time_delay( dec, ra, ifos[1], ifos[0], coord=opts.coord, tgeocent=opts.gps, degrees=False )
        ylim = ax.get_ylim()
        ax.plot( [dt*1e3]*2, ylim, color=opts.time_delay_color, alpha=opts.time_delay_alpha )
        ax.set_ylim(ylim)

    ### decorate axis
    ax.set_xlim(xmin=xmin*1e3, xmax=xmax*1e3)

    aX = ax.twiny()
    aX.set_xticklabels(["$%.1f^\circ$"%(np.arccos(tick*1e-3/maxDt)*180/np.pi) for tick in ax.get_xticks()])
    aX.set_xlim(ax.get_xlim())

    ax.set_xlabel(r'$\Delta t_{%s}\ [\mathrm{ms}]$'%(ifos))
    ax.set_ylabel(r'$p(\Delta t_{%s}|\mathrm{data})$'%(ifos))

    if opts.no_yticks:
        ax.set_yticklabels([])

    ax.grid(True, which='both')

#    ax.legend(loc='best')

    for figtype in opts.figtype:
        figname = "%s/dT_%s%s.%s"%(opts.output_dir, ifos, opts.tag, figtype)
        if opts.verbose:
            print figname
        fig.savefig(figname, dpi=opts.dpi)
    plt.close(fig)
