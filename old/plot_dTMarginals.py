#!/usr/bin/python
usage       = "plot_dTMarginals.py [--options] label,FITS label,FITS ..."
description = "plot the marigals for time of flight between the two specified IFOs"
author      = "R. Essick (reed.essick@ligo.org)"

#-------------------------------------------------

import stats
import triangulate 

import healpy as hp
import numpy as np

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
plt.rcParams.update({'font.family':'serif'})

from optparse import OptionParser

#-------------------------------------------------

twopi = 2*np.pi
pi2 = 0.5*np.pi

colors = ['b', 'r', 'g', 'm', 'c', 'k', 'y']

axpos = [0.10, 0.10, 0.85, 0.85]

#-------------------------------------------------

def gen_fig_ax( figind, figwidth=9, figheight=5 ):
    '''
    generates figure and axis in the set-up we prefer
    '''
    fig = plt.figure(figind, figsize=(figwidth, figheight) )
    ax = fig.add_axes(axpos)
    ax.grid( True )

    return fig, ax

def post2marg( post, ifos, sampDt, coord='C', gps=None ):
    '''
    compute the marginal distribution for time-of-flight from the full posterior over the sky
    '''
    nside = hp.npix2nside(len(post))
    dAng = hp.nside2pixarea(nside)**0.5

    maxDt = sampDt[-1]

    y, x = hp.pix2ang(nside, np.arange(len(post)))
    if coord=="C": ### convert theta -> dec
        y = pi2 - y

    ### compute time-of-flight
    dt = triangulate.time_delay( y, x, ifos[1], ifos[0], coord=coord, tgeocent=gps, degrees=False )

    ### build up KDE
    sigma = maxDt * ( 1 - dt/maxDt )**0.5 * dAng

    spacer = np.ones_like( sampDt )

    big_post   = np.outer( spacer, post )
    big_dt     = np.outer( spacer, dt )
    big_sigma  = np.outer( spacer, sigma )
    big_sampDt = np.outer( sampDt, np.ones_like(dt) )

    kde = np.sum( big_post * np.exp( -( (big_dt - big_sampDt)/big_sigma )**2 ) / (twopi**0.5 * big_sigma), axis=-1 )
    kde /= np.sum(kde)

    return kde

def plot( ax, sampDt, marg, color='k', label=None, xlim_dB=-20 ):
    '''
    plot data
    '''
    ax.plot( sampDt*1e3, marg*1e-3, label=label, color=color )

    subset = sampDt[marg > 10**(xlim_dB/10)*np.max(marg)]
    ax.set_xlim( xmin=min(xmin, np.min(subset)), xmax=max(xmax, np.max(subset)) )

def set_xlim( ax, xmin=None, xmax=None ):
    '''
    sets the xlimits
    '''
    if xmin!=None:
        ax.set_xlim( xmin=xmin )
    if xmax!=None:
        ax.set_xlim( xmax=xmax )

def annotate( ax, SRCs, IFOs, maxDt, color='k', alpha=1.0, coord='C', gps=None, degrees=False ):
    '''
    annotate plot
    '''
    for dec, ra in SRCs:
        if degrees:
            dec *= triangulate.deg2rad
            ra  *= triangulate.deg2rad

        dt = triangulate.time_delay( dec, ra, ifos[1], ifos[0], coord=coord, tgeocent=gps, degrees=False )
        ylim = ax.get_ylim()
        ax.plot( [dt*1e3]*2, ylim, color=color, alpha=alpha )
        ax.set_ylim(ylim)

    aX = ax.twiny()
    aX.set_xticklabels(["$%.1f^\circ$"%(np.arccos(tick*1e-3/maxDt)*180/np.pi) for tick in ax.get_xticks()])
    aX.set_xlim(ax.get_xlim())

#------------------------

def gen_sampDt( IFOs, Nsamp=501 ):
    '''
    generate sample points for KDE generation
    '''
    t, p = triangulate.line_of_sight( ifos[0], ifos[1], coord='E' )
    maxDt = np.abs( triangulate.time_delay( t, p, ifos[0], ifos[1], coord='E' ) )
    return np.linspace(-maxDt, maxDt, Nsamp)

#-------------------------------------------------

if __name__=="__main__":

    ### parse arguments

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

    #---------------------------------------------

    if opts.verbose:
        print "reading maps"
    maps = {}
    for arg in args:
        label, fits = arg.split(",")
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
        sampDt = gen_sampDt( ifos, Nsamp=opts.Nsamp )
        maxDt = sampDt[-1]

        ### compute and plot marginals
        fig, ax = gen_fig_ax( figind, figwidth=opts.figwidth, figheight=opts.figheight )
        figind += 1

        ax.set_xlim(xmin=maxDt, xmax=-maxDt)

        for cind, label in enumerate(labels):
            if opts.verbose:
                print "    ", label

            post  = maps[label]['post']

            kde = post2marg( post, ifos, sampDt, coord=opts.coord, gps=opts.gps ) 

            ### plot
            color = colors[cind%len(colors)]
            plot( sampDT, kde, label=label, color=color, xlim_dB=opts.xlim_dB )
            fig.text(0.10+0.05, 0.95-0.05*(figind-1), label, color=color, ha='left', va='top')

        ### annotate the plot
        annotate( ax, 
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

        ax.legend(loc='best')

        if opts.no_yticks:
            ax.set_yticklabels([])

        if opts.transparent:
            stack_fig.patch.set_alpha(0.)
            stack_ax.patch.set_alpha(0.)
            stack_ax.set_alpha(0.)

        ### save figure
        for figtype in opts.figtype:
            figname = "%s/dT_%s%s.%s"%(opts.output_dir, ifos, opts.tag, figtype)
            if opts.verbose:
                print figname
            fig.savefig(figname, dpi=opts.dpi)
        plt.close(fig)
