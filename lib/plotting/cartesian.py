description = "a module that houses convenient functions for plotting skymaps and skymap-related things, focusing on mollweide projections"
author      = "reed.essick@ligo.org"

#-------------------------------------------------

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
plt.rcParams.update({'font.family':'serif', 'text.usetex':True})

import stats
import triangulate

import healpy as hp
import numpy as np

#-------------------------------------------------

twopi = 2*np.pi
pi2 = 0.5*np.pi

dT_axpos = [0.10, 0.10, 0.85, 0.85]

primetime  = [0.10, 0.10, 0.65, 0.55]
right_proj = [0.76, 0.10, 0.19, 0.55]
top_proj   = [0.10, 0.66, 0.65, 0.29]

#-------------------------------------------------

### actual plotting and figure manipulation

def genDT_fig_ax( figind, figwidth=9, figheight=5 ):
    '''
    generates figure and axis in the set-up we prefer
    '''
    fig = plt.figure(figind, figsize=(figwidth, figheight) )
    ax = fig.add_axes(dT_axpos)
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
    xmin, xmax = ax.get_xlim()

    ax.plot( sampDt*1e3, marg*1e-3, label=label, color=color )

    subset = sampDt[marg > 10**(xlim_dB/10)*np.max(marg)]*1e3
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

#-----------

def genHist_fig_ax( figind, figwidth=9, figheight=5 ):
    fig = plt.figure(figind, figsize=(figwidth, figheight) )

    ax    = fig.add_axes(primetime)
    rproj = fig.add_axes(right_proj)
    tproj = fig.add_axes(top_proj)

    ax.grid(True)
    rproj.grid(True)
    tproj.grid(True)

    return fig, ax, rproj, tproj

def histogram2d( theta, phi, ax, rproj, tproj, weights=None, Nbins=250, color='b', log=False, contour=False, cmap="jet" ):
    """
    custom plot for skymap analysis
    """
    ### define binning
    theta_bins = np.linspace(0, np.pi, Nbins+1)
    theta_dots = 0.5*(theta_bins[:-1]+theta_bins[1:]) * 180 / np.pi

    phi_bins = np.linspace(-np.pi, np.pi, Nbins+1)
    phi_dots = 0.5*(phi_bins[:-1]+phi_bins[1:]) * 180 / np.pi

    ### 2D histogram
    tp_count = np.histogram2d( theta, phi, bins=(theta_bins, phi_bins), weights=weights )[0]
    if contour:
        if log:
            ax.contour( phi_dots, theta_dots, np.log10(tp_count), colors=color, alpha=0.5 )
        else:
            ax.contour( phi_dots, theta_dots, tp_count, colors=color, alpha=0.5 )
    else:
        im = matplotlib.image.NonUniformImage( ax, interpolation='bilinear', cmap=plt.get_cmap(cmap) )
        if log:
            im.set_data( phi_dots, theta_dots[::-1], np.log10(tp_count[::-1,:]) )
        else:
            im.set_data( phi_dots, theta_dots[::-1], tp_count[::-1,:] )
        im.set_alpha( 0.001 )
        ax.images.append( im )

    ax.set_xlim( xmin=-180.0, xmax=180.0 )
    ax.set_ylim( ymin=180.0, ymax=0.0 )
    ax.set_xlabel( "$\phi$" )
    ax.set_ylabel( "$\\theta$" )

    ax.set_xticks( np.arange(-180, 180, 10), minor=True )
    ax.set_yticks( np.arange(0, 180, 5), minor=True )

    ### theta projection
    theta_count = np.histogram( theta, bins=theta_bins, weights=weights )[0]
    rproj.plot( theta_count, theta_dots[::-1], color=color )
    rproj.set_ylim( ymin=0.0, ymax=180.0 )
    rproj.set_xlabel( "$p(\\theta)$" )
    plt.setp(rproj.get_yticklabels(), visible=False)
    rproj.set_yticks( np.arange(0, 180, 5), minor=True )
    if log:
        rproj.set_xscale('log')

    ### phi proejection
    phi_count = np.histogram( phi, bins=phi_bins, weights=weights )[0]
    tproj.plot( phi_dots, phi_count, color=color )
    tproj.set_xlim( xmin=-180.0, xmax=180.0 )
    tproj.set_ylabel( "$p(\phi)$" )
    plt.setp(tproj.get_xticklabels(), visible=False)
    tproj.set_xticks( np.arange(-180, 180, 10), minor=True )
    if log:
        tproj.set_yscale('log')

#------------------------

### data preparation

def gen_sampDt( IFOs, Nsamp=501 ):
    '''
    generate sample points for KDE generation
    '''
    t, p = triangulate.line_of_sight( IFOs[0], IFOs[1], coord='E' )
    maxDt = np.abs( triangulate.time_delay( t, p, IFOs[0], IFOs[1], coord='E' ) )
    return np.linspace(-maxDt, maxDt, Nsamp)
