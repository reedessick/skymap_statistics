description = "a module that houses convenient functions for plotting skymaps and skymap-related things, focusing on mollweide projections"
author      = "reed.essick@ligo.org"

#-------------------------------------------------

import os
import json

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
try:
        from lalinference import plot as lalinf_plot
        from lalinference import cmap
except:
        raise StandardError("Could not import lalinference.plot")
plt.rcParams.update({'font.family':'serif', 'text.usetex':True })

import numpy as np
import healpy as hp

import triangulate

#-------------------------------------------------

axpos = [0.03, 0.03, 0.94, 0.94]

twopi = 2*np.pi

#-------------------------------------------------

### actual plotting and figure manipulation

def gen_fig_ax( figind, figheight=5, figwidth=9, projection=None ):
    '''
    generates figure and axis in the set-up we prefer
    '''
    fig = plt.figure(figind, figsize=(figwidth, figheight) )
    if projection:
        ax = fig.add_axes(axpos, projection=projection)
    else:
        ax = fig.add_axes(axpos)
    ax.grid( True )

    return fig, ax

def annotate( ax, projection=None, line_of_sight=[], line_of_sight_color='k', zenith=[], zenith_color='k', time_delay=[], time_delay_color='k', time_delay_alpha=1.0, marker_Dec_RA=[], marker='o', marker_color='k', marker_size=4, marker_edgewidth=1, marker_alpha=1.0, continents=False, continents_color='k', continents_alpha=1.0 ):
    '''
    annotates the mollweide projection
    '''
    ### plot line-of-sight markers and text
    for ifos, (y,x), (Y,X) in line_of_sight:
        if projection=="mollweide":
            if x > np.pi:
                x -= twopi
            if X > np.pi:
                X -= twopi

        ax.plot( x, y, color=line_of_sight_color, markeredgecolor=line_of_sight_color, marker='o', markersize=2 )
        ax.plot( X, Y, color=line_of_sight_color, markeredgecolor=line_of_sight_color, marker='o', markersize=2 )

        ax.text( x, y, " %s-%s"%(ifos[1],ifos[0]), ha='left', va='bottom', color=line_of_sight_color )
        ax.text( X, Y, " %s-%s"%(ifos[0],ifos[1]), ha='left', va='bottom', color=line_of_sight_color )

    ### plot zenith markers
    for ifo, (y,x), (Y,X) in zenith:
        if projection=="mollweide":
            if x > np.pi:
                x -= twopi
            if X > np.pi:
                X -= twopi

        ax.plot( x, y, color=zenith_color, markeredgecolor=zenith_color, marker='s', markersize=2 )
        ax.plot( X, Y, color=zenith_color, markeredgecolor=zenith_color, marker='s', markersize=2 )

        ax.text( x, y, " "+ifo+"+", ha='left', va='bottom', color=zenith_color )
        ax.text( X, Y, " "+ifo+"-", ha='left', va='bottom', color=zenith_color )

    ### plot time-delay loci
    for y, x in time_delay:
        ax.plot( x, y, color=time_delay_color, alpha=time_delay_alpha )

    ### plot general markers
    for dec, ra in marker_Dec_RA:
        ax.plot( ra, dec,
                 linestyle='none',
                 marker=marker,
                 markerfacecolor='none',
                 markeredgecolor=marker_color,
                 markersize=marker_size,
                 markeredgewidth=marker_edge_width,
                 alpha=marker_alpha )

    ### add continents
    if continents:
        geojson_filename = os.path.join(os.path.dirname(lalinf_plot.__file__), 'ne_simplified_coastline.json')
        file_obj = open(geojson_filename, 'r')
        geojson = json.load(file_obj)
        file_obj.close()

        for shape in geojson['geometries']:
            verts = np.deg2rad(shape['coordinates'])
            ax.plot( verts[:, 0], verts[:, 1], color=continents_color, linewidth=0.5, alpha=continents_alpha )

def heatmap( post, ax, color_map='OrRd' ):
    '''
    generate mollweide projection of heatmap with requested annotations
    '''
    plt.sca( ax )
    lalinf_plot.healpix_heatmap( post, cmap=plt.get_cmap(color_map) )

def contour( post, ax, levels=[0.1, 0.5, 0.9], alpha=1.0, colors='b', linewidths=1 ):
    '''
    generate mollweide projection of contours with requested annotations
    '''
    cpost = np.empty(post.shape)
    indecies = np.argsort(post)[::-1]
    cpost[indecies] = np.cumsum(post[indecies])

    plt.sca( ax )
    lalinf_plot.healpix_contour( cpost,
                                 levels=levels,
                                 alpha=alpha,
                                 colors=colors,
                                 linewidths=linewidths )

#-------------------------------------------------

### data preparation

def gen_line_of_sight( IFOs, coord='C', gps=None ):
    '''
    computes line-of-sight directions and returns them in a friendly format for plotting
    '''
    line_of_sight = []
    for ifos in IFOs:
        y, x = triangulate.line_of_sight(ifos[1], ifos[0], coord=coord, tgeocent=gps, degrees=False)
        X, Y = triangulate.antipode( x, y, coord=coord, degrees=False)
        if coord=="E": ### convert theta->dec
            y = 0.5*np.pi - y
            Y = 0.5*np.pi - Y
        line_of_sight.append( (ifos, (y,x), (Y,X)) )
    return line_of_sight

def gen_zenith( IFOs, coord='C', gps=None ):
    '''
    computes zenith directions and returns them in a friendly format for plotting
    '''
    zenith = []
    for ifo in IFOs:
        y, x = triangulate.overhead(ifo, coord=coord, tgeocent=gps, degrees=False)
        X, Y = triangulate.antipode( x, y, coord=coord, degrees=False)
        if coord=="E": ### convert theta->dec
            y = 0.5*np.pi - y
            Y = 0.5*np.pi - Y
        zenith.append( (ifo, (y,x), (Y,X)) )
    return zenith

def gen_time_delay( SRCs, IFOs, coord='C', gps=None, degrees=False ):
    '''
    computes locus of points with constant time-delay between IFOs and returns them in a friendly format for plotting
    '''
    time_delay = []
    for dec, ra in SRCs:
        if degrees:
            dec *= triangulate.deg2rad
            ra *= triangulate.deg2rad

        for ifos in IFOs:
            dt = triangulate.time_delay( dec, ra, ifos[1], ifos[0], coord=coord, tgeocent=gps, degrees=False )
            y, x = triangulate.time_delay_locus( dt, ifos[1], ifos[0], coord=coord, tgeocent=gps, degrees=False )

            if coord=="E": ### convert theta-> dec
                y = 0.5*np.pi - y

                x[x>np.pi] -= 2*np.pi ### ensure that everything is between -pi and pi

            ### find big jumps in azimuthal angle and split up the plotting jobs
            d = np.concatenate( ([0],np.nonzero(np.abs(x[1:]-x[:-1])>np.pi)[0]+1,[len(x)]) )
            for istart, iend in zip(d[:-1], d[1:]):
                time_delay.append( (y[istart:iend], x[istart:iend]) )

    return time_delay

def gen_marker_Dec_RA( SRCs, coord='C', gps=None, degrees=False ):
    '''
    returns marker positions in a friendly format for plotting
    '''
    marker_Dec_RA = []
    if coord=="E":
        for dec, ra in SRCs:
            if degrees:
                dec *= triangulate.deg2rad
                ra *= triangulate.deg2rad
            if ra > np.pi:
                ra -= 2*np.pi
            marker_Dec_RA.append( (0.5*np.pi-dec, ra) )
    else: ### coord=="C"
        for dec, ra in SRCs:
            if degrees:
                dec *= triangulate.deg2rad
                ra *= triangulate.deg2rad
#            if ra > np.pi:
#                ra -= 2*np.pi
            marker_Dec_RA.append( (dec, ra) )

    return marker_Dec_RA
