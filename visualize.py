description = "module to hold plotting functions"
author = "Reed Essick (reed.essick@ligo.org)"

import numpy as np

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
plt.rcParams.update( {"text.usetex":True} )

#=================================================

def histogram2d( theta, phi, weights=None, figax=None, Nbins=500, color='b', log=False ):
    """
    custom plot for skymap analysis
    """
    primetime = [0.10, 0.10, 0.65, 0.55]
    right_proj = [0.76, 0.10, 0.19, 0.55]
    top_proj = [0.10, 0.66, 0.65, 0.29]

    if figax:
        fig, ax = figax
    else:
#        fig = plt.figure(projection="astro mollweide")
        fig = plt.figure( figsize=(10,8) )
        ax = [fig.add_axes(primetime), fig.add_axes(right_proj), fig.add_axes(top_proj)]

    pt, rp, tp = ax

    ### define binning
    theta_bins = np.linspace(0, np.pi, Nbins+1)
    theta_dots = 0.5*(theta_bins[:-1]+theta_bins[1:]) / np.pi

    phi_bins = np.linspace(-np.pi, np.pi, Nbins+1)
    phi_dots = 0.5*(phi_bins[:-1]+phi_bins[1:]) / np.pi

    ### 2D histogram
    tp_count = np.histogram2d( theta, phi, bins=(theta_bins, phi_bins), weights=weights )[0]

    im = matplotlib.image.NonUniformImage( pt, interpolation='bilinear')
    if log:
        im.set_data( phi_dots, theta_dots[::-1], np.log10(tp_count[::-1,:]) )
    else:
        im.set_data( phi_dots, theta_dots[::-1], tp_count[::-1,:] )
    im.set_alpha( 0.001 )
    pt.images.append( im )
    pt.set_xlim( xmin=-1.0, xmax=1.0 )
    pt.set_ylim( ymin=1.0, ymax=0.0 )
    pt.set_xlabel( "$\phi$" )
    pt.set_ylabel( "$\\theta$" )

    ### theta projection
    theta_count = np.histogram( theta, bins=theta_bins, weights=weights )[0]
    rp.plot( theta_count, theta_dots[::-1], color=color )
    rp.set_ylim( ymin=0.0, ymax=1.0 )
    rp.set_xlabel( "$p(\\theta)$" )
    plt.setp(rp.get_yticklabels(), visible=False)
    if log:
        rp.set_xscale('log')

    ### phi proejection
    phi_count = np.histogram( phi, bins=phi_bins, weights=weights )[0]
    tp.plot( phi_dots, phi_count, color=color )
    tp.set_xlim( xmin=-1.0, xmax=1.0 )
    tp.set_ylabel( "$p(\phi)$" )
    plt.setp(tp.get_xticklabels(), visible=False)
    if log:
        tp.set_yscale('log')

    return fig, ax

def corner( theta, phi, weights=None ):
    """
    use the corner module
    """
    import corner

    fig = corner.corner( np.array([theta, phi]).transpose(), 
                         weights=weights, 
                         show_titles=True,
                         truths=None,
                         quantiles = [0.1, 0.5, 0.9],
                         labels = ["$\\theta$", "$\phi$"]
                       )

    return fig, None
