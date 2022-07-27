description = "a module that houses convenient functions for plotting skymaps and skymap-related things, focusing on cartesian projections"
author      = "reed.essick@ligo.org"

#-------------------------------------------------

import os
import json

import healpy as hp
import numpy as np

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
plt.rcParams.update({'font.family':'serif', 'text.usetex':True})

### non-standard libraries
from skymap_statistics import stats
from skymap_statistics import triangulate

#-------------------------------------------------

twopi = 2*np.pi
pi2 = 0.5*np.pi
deg2rad = np.pi/180
rad2deg = 1./deg2rad

dT_axpos = [0.10, 0.10, 0.85, 0.85]

primetime  = [0.10, 0.10, 0.65, 0.55]
right_proj = [0.76, 0.10, 0.19, 0.55]
top_proj   = [0.10, 0.66, 0.65, 0.29]

cr_axpos = [0.10, 0.10, 0.85, 0.85]

#-------------------------------------------------

### actual plotting and figure manipulation

def gen_fig_ax( figind, figwidth=6, figheight=5 ):
    '''
    generates figure and axis in the set-up we prefer
    '''
    return genCR_fig_ax(figind, figwidth=figwidth, figheigh=figheight)

def genCR_fig_ax( figind, figwidth=6, figheight=5.5, grid=True ):
    '''
    generates figure and axis in the set-up we prefer
    '''
    fig = plt.figure(figind, figsize=(figwidth, figheight))
    ax = fig.add_axes(cr_axpos)
    ax.grid(grid, which='both')

    return fig, ax

def genDT_fig_ax( figind, figwidth=9, figheight=5, grid=True ):
    '''
    generates figure and axis in the set-up we prefer
    '''
    fig = plt.figure(figind, figsize=(figwidth, figheight) )
    ax = fig.add_axes(dT_axpos)
    ax.grid(grid, which='both')

    return fig, ax

def genHist_fig_ax( figind, figwidth=9, figheight=5, grid=True ):
    fig = plt.figure(figind, figsize=(figwidth, figheight) )

    ax    = fig.add_axes(primetime)
    rproj = fig.add_axes(right_proj)
    tproj = fig.add_axes(top_proj)

    ax.grid(grid, which='both')
    rproj.grid(grid, which='both')
    tproj.grid(grid, which='both')

    return fig, ax, rproj, tproj

def gen_limits(minX, maxX, minY, maxY, coord='C', degrees=False):
    '''
    set up the xlimits and ylimits for plotting
    return xlim, ylim
    '''
    ### insure we have everything defined
    if minX==None:
        minX = 0

    if maxX==None:
        maxX = 360 if degrees else twopi

    if minY==None:
        if coord=='C':
            minY = -90 if degrees else -pi2
        else:
            minY = 0

    if maxY==None:
        if coord=='C':
            maxY = 90 if degrees else pi2
        else:
            maxY = 180 if degrees else np.pi

    ### insure everyting is in radians
    if degrees:
        minX *= deg2rad
        maxX *= deg2rad
        minY *= deg2rad
        maxY *= deg2rad

    ### insure we have ordering correct based on coordinate system
    if coord=='C': ### RA increases to the left, 
        minX, maxX = maxX, minX ### RA increases to the left...

    elif coord=='E':
        minY, maxY = pi2 - np.array([maxY, minY]) ### change theta -> Dec
        
    else:
        raise ValueError( 'coord=%s not understood'%coord )

    ### insure evertyhing wraps correctly
    if minX!=twopi:
        minX = minX%(twopi)
    if maxX!=twopi:
        maxX = maxX%(twopi)
    assert -pi2<=minY<=pi2, 'minY out of bounds!'
    assert -pi2<=maxY<=pi2, 'maxY out of bounds!'

    return (minX, maxX), (minY, maxY)

#------------------------

def heatmap( post, ax, xlim, ylim, color_map='jet', Npts=1001, colorbar=False, colorbar_label='' ):
    '''
    take a HEALPix posterior and make a cartesian heatmap
    '''
    cart = post2cart(post, xlim, ylim, Npts=Npts)
    ax.imshow(cart, 
        interpolation='bilinear', 
        origin='lower', 
        extent=(xlim[0], xlim[1], ylim[0], ylim[1]),
        aspect='auto',
        cmap=plt.get_cmap(color_map),
    )
    if colorbar:
        cb = plt.colorbar(ax=ax)
        cb.set_label(colorbar_label)

def contour( post, ax, xlim, ylim, levels=[0.1, 0.5, 0.9], alpha=1.0, colors='b', linewidths=1, Npts=1001 ):
    '''
    take a HEALPix posterior and make a cartesian contour plot
    '''
    cpost = np.empty(post.shape)
    indecies = np.argsort(post)[::-1]
    cpost[indecies] = np.cumsum(post[indecies])

    cart = post2cart(cpost, xlim, ylim, Npts=Npts)
    ax.contour(cart, 
        levels, 
        alpha=alpha, 
        colors=colors, 
        linewidths=linewidths,
        origin='lower',
        extent=(xlim[0], xlim[1], ylim[0], ylim[1]),
    )

def plot_dT( ax, sampDt, marg, color='k', label=None, xlim_dB=-20 ):
    '''
    plot data
    '''
    xmin, xmax = ax.get_xlim()

    ax.plot( sampDt*1e3, marg*1e-3, label=label, color=color )

    subset = sampDt[marg > 10**(xlim_dB/10)*np.max(marg)]*1e3
    ax.set_xlim( xmin=min(xmin, np.min(subset)), xmax=max(xmax, np.max(subset)) )

def histogram2d( theta, phi, ax, rproj, tproj, weights=None, Nbins=250, color='b', log=False, contour=False, levels=[0.1, 0.5, 0.9], cmap="jet", alpha=1.0 ):
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
        ### convert to cumulative count for countours!
        shape = tp_count.shape

        tp_count = tp_count.flatten()
        tp_cum = np.empty(tp_count.shape)
        order = np.argsort(tp_count)[::-1]
        tp_cum[order] = np.cumsum( tp_count[order] )

        tp_cum = np.reshape(tp_cum, shape)

        if log:
            ax.contour( phi_dots, theta_dots, np.log10(tp_cum), colors=color, alpha=alpha, levels=levels )
        else:
            ax.contour( phi_dots, theta_dots, tp_cum, colors=color, alpha=alpha, levels=levels )
    else:
        im = matplotlib.image.NonUniformImage( ax, interpolation='bilinear', cmap=plt.get_cmap(cmap) )
        if log:
            im.set_data( phi_dots, theta_dots[::-1], np.log10(tp_count[::-1,:]) )
        else:
            im.set_data( phi_dots, theta_dots[::-1], tp_count[::-1,:] )
        im.set_alpha( alpha )
        ax.images.append( im )

    ax.set_xlim( xmin=-180.0, xmax=180.0 )
    ax.set_ylim( ymin=180.0, ymax=0.0 )
    ax.set_xlabel( "$\phi$" )
    ax.set_ylabel( "$\\theta$" )

    ax.set_xticks( np.arange(-180, 180, 10), minor=True )
    ax.set_yticks( np.arange(0, 180, 5), minor=True )

    ### theta projection
    theta_count = np.histogram( theta, bins=theta_bins, weights=weights )[0]
    rproj.plot( theta_count, theta_dots[::-1], color=color, alpha=alpha)
    rproj.set_ylim( ymin=0.0, ymax=180.0 )
    rproj.set_xlabel( "$p(\\theta)$" )
    plt.setp(rproj.get_yticklabels(), visible=False)
    rproj.set_yticks( np.arange(0, 180, 5), minor=True )
    if log:
        rproj.set_xscale('log')

    ### phi proejection
    phi_count = np.histogram( phi, bins=phi_bins, weights=weights )[0]
    tproj.plot( phi_dots, phi_count, color=color, alpha=alpha )
    tproj.set_xlim( xmin=-180.0, xmax=180.0 )
    tproj.set_ylabel( "$p(\phi)$" )
    plt.setp(tproj.get_xticklabels(), visible=False)
    tproj.set_xticks( np.arange(-180, 180, 10), minor=True )
    if log:
        tproj.set_yscale('log')

#------------------------

def set_xlim( ax, xmin=None, xmax=None ):
    '''
    sets the xlimits
    '''
    if xmin!=None:
        ax.set_xlim( xmin=xmin )
    if xmax!=None:
        ax.set_xlim( xmax=xmax )

def set_ylim( ax, ymin=None, ymax=None ):
    '''
    set the ylimits
    '''
    if ymin!=None:
        ax.set_ylim( ymin=ymin )
    if ymax!=None:
        ax.set_ylim( ymax=ymax )

def set_lim( ax, xmin=None, xmax=None, ymin=None, ymax=None ):
    '''
    set all limits
    '''
    set_xlim(ax, xmin=xmin, xmax=xmax)
    set_ylim(ax, ymin=ymin, ymax=ymax)

def set_labels( ax, coord='C' ):
    '''
    set the axis labels (and modify ticklabels) according to which coordinate system is specified
    '''
    ax.set_xticklabels(['$%.1f^\circ$'%(_*rad2deg) for _ in ax.get_xticks()])

    ### set labels
    if coord=="C":
        ax.set_xlabel('RA')
        ax.set_ylabel('Dec')
        ax.set_yticklabels(['$%.1f^\circ$'%(_*rad2deg) for _ in ax.get_yticks()])

    elif coord=="E":
        ax.set_xlabel('$\phi$')
        ax.set_ylabel(r'$\theta$')
        ax.set_yticklabels(['$%.1f^\circ$'%((pi2-_)*rad2deg) for _ in ax.get_yticks()])

    else:
        raise ValueError( 'coord=%s not understood'%coord )

def annotate( ax, line_of_sight=[], line_of_sight_color='k', line_of_sight_fontsize=8, line_of_sight_markersize=2, zenith=[], zenith_color='k', zenith_fontsize=8, zenith_markersize=2, time_delay=[], time_delay_color='k', time_delay_alpha=1.0, time_delay_linestyle='solid', marker_Dec_RA=[], marker='o', marker_color='k', marker_size=4, marker_edgewidth=1, marker_alpha=1.0, continents=[], continents_color='k', continents_alpha=1.0, constellations=[], constellations_color='k', constellations_alpha=1.0, stars=[], stars_color='k', stars_alpha=1.0, arms=[], arms_color='k', arms_linewidth=1, arms_alpha=1.0, constellation_boundaries=[], constellation_boundaries_color='k', constellation_boundaries_alpha=1.0, constellation_centers=[], constellation_centers_color='k', constellation_centers_alpha=1.0, constellation_centers_fontsize=8 ):
    '''
    annotate cartesian projection
    '''
    ### plot line-of-sight markers and text
    for ifos, (y,x), (Y,X) in line_of_sight:
        ax.plot( x, y, color=line_of_sight_color, markeredgecolor=line_of_sight_color, marker='o', markersize=line_of_sight_markersize )
        ax.plot( X, Y, color=line_of_sight_color, markeredgecolor=line_of_sight_color, marker='o', markersize=line_of_sight_markersize )

        ax.text( x, y, " %s-%s"%(ifos[1],ifos[0]), ha='left', va='bottom', color=line_of_sight_color, fontsize=line_of_sight_fontsize )
        ax.text( X, Y, " %s-%s"%(ifos[0],ifos[1]), ha='left', va='bottom', color=line_of_sight_color, fontsize=line_of_sight_fontsize )

    ### plot zenith markers
    for ifo, (y,x), (Y,X) in zenith:
        ax.plot( x, y, color=zenith_color, markeredgecolor=zenith_color, marker='s', markersize=zenith_markersize )
        ax.plot( X, Y, color=zenith_color, markeredgecolor=zenith_color, marker='s', markersize=zenith_markersize )

        ax.text( x, y, " "+ifo+"+", ha='left', va='bottom', color=zenith_color, fontsize=zenith_fontsize )
        ax.text( X, Y, " "+ifo+"-", ha='left', va='bottom', color=zenith_color, fontsize=zenith_fontsize )

    ### plot time-delay loci
    for y, x in time_delay:
        ax.plot( x, y, color=time_delay_color, alpha=time_delay_alpha, linestyle=time_delay_linestyle )
        ax.plot( x+twopi, y, color=time_delay_color, alpha=time_delay_alpha, linestyle=time_delay_linestyle )

    ### plot general markers
    for dec, ra in marker_Dec_RA:
        ax.plot( ra, dec,
                 linestyle='none',
                 marker=marker,
                 markerfacecolor='none',
                 markeredgecolor=marker_color,
                 markersize=marker_size,
                 markeredgewidth=marker_edgewidth,
                 alpha=marker_alpha )

    ### add continents
    for verts in continents: ### plot repeatedly to account for periodicity
        ax.plot( verts[:, 0], verts[:, 1], color=continents_color, linewidth=0.5, alpha=continents_alpha )
        ax.plot( verts[:, 0]+twopi, verts[:, 1], color=continents_color, linewidth=0.5, alpha=continents_alpha )
        ax.plot( verts[:, 0]-twopi, verts[:, 1], color=continents_color, linewidth=0.5, alpha=continents_alpha )

    ### add arms
    for x, y in arms:
        ax.plot( x, y, color=arms_color, linewidth=arms_linewidth, alpha=arms_alpha )
        ax.plot( x+twopi, y, color=arms_color, linewidth=arms_linewidth, alpha=arms_alpha )
        ax.plot( x-twopi, y, color=arms_color, linewidth=arms_linewidth, alpha=arms_alpha )

    ### add constellations
    for shape in constellations:
        ax.plot( shape[:,0], shape[:,1], color=constellations_color, linewidth=0.5, alpha=constellations_alpha )
        ax.plot( shape[:,0]+twopi, shape[:,1], color=constellations_color, linewidth=0.5, alpha=constellations_alpha )
        ax.plot( shape[:,0]-twopi, shape[:,1], color=constellations_color, linewidth=0.5, alpha=constellations_alpha )

    ### add stars
    for x, y, mag in stars:
        markersize=max(1, 5-mag) ### FIXME: hard coded...bad?
        if x < 0: 
            x += twopi
        elif x > twopi: 
            x -= twopi
        ax.plot(x, y, markersize=markersize, marker='o', markerfacecolor=stars_color, markeredgecolor='none', alpha=stars_alpha)

    ### add constellation boundaries
    for x, y in constellation_boundaries:
        ax.plot(x, y, linewidth=0.5, color=constellation_boundaries_color, alpha=constellation_boundaries_alpha)
        ax.plot(x-twopi, y, linewidth=0.5, color=constellation_boundaries_color, alpha=constellation_boundaries_alpha)
        ax.plot(x+twopi, y, linewidth=0.5, color=constellation_boundaries_color, alpha=constellation_boundaries_alpha)

    ### add constellation centers
    for x, y, name in constellation_centers:
        if x < 0:
            x += twopi
        elif x > twopi:
            x -= twopi
        ax.text(
            x, y, name, 
            ha='center', 
            va='center', 
            alpha=constellation_centers_alpha, 
            color=constellation_centers_color, 
            fontsize=constellation_centers_fontsize,
        )

def annotateDT( ax, SRCs=[], IFOs='HL', color='k', alpha=1.0, coord='C', gps=None, degrees=False, twiny=True ):
    '''
    annotate time-delay plot
    '''
    for dec, ra in SRCs:
        if degrees:
            dec *= triangulate.deg2rad
            ra  *= triangulate.deg2rad

        dt = triangulate.time_delay( dec, ra, IFOs[1], IFOs[0], coord=coord, tgeocent=gps, degrees=False )
        ylim = ax.get_ylim()
        ax.plot( [dt*1e3]*2, ylim, color=color, alpha=alpha, linestyle='--' )
        ax.set_ylim(ylim)

    if twiny:
        t, p = triangulate.line_of_sight( IFOs[0], IFOs[1], coord='E' )
        maxDt = np.abs( triangulate.time_delay( t, p, IFOs[0], IFOs[1], coord='E' ) )
       
        aX = ax.twiny()
        aX.set_xticklabels(["$%.1f^\circ$"%(np.arccos(tick*1e-3/maxDt)*180/np.pi) for tick in ax.get_xticks()])
        aX.set_xlim(ax.get_xlim())

#------------------------

### data preparation

def gen_sampDt( IFOs, Nsamp=501 ):
    '''
    generate sample points for KDE generation
    '''
    t, p = triangulate.line_of_sight( IFOs[0], IFOs[1], coord='E' )
    maxDt = np.abs( triangulate.time_delay( t, p, IFOs[0], IFOs[1], coord='E' ) )
    return np.linspace(-maxDt, maxDt, Nsamp)

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
    
    kde = np.array([ np.sum(post * np.exp( -( (dt - k)/sigma )**2 ) / (twopi**0.5 * sigma)) for k in sampDt])
    kde /= np.sum(kde)
    return kde

def post2cart( post, ra_lim, dec_lim, Npts=1001 ):
    '''
    convert a healpix posterior into a cartesian mapping
    '''
    #                                                         convert dec -> theta
    PHI, THETA = np.meshgrid(np.linspace(*ra_lim, num=Npts), pi2 - np.linspace(*dec_lim, num=Npts))
    PHI = PHI.flatten()
    THETA = THETA.flatten()
    return post[hp.ang2pix(hp.npix2nside(len(post)), THETA, PHI)].reshape((Npts, Npts))
