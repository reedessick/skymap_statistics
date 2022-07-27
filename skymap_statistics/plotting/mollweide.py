description = "a module that houses convenient functions for plotting skymaps and skymap-related things, focusing on mollweide projections"
author      = "reed.essick@ligo.org"

#-------------------------------------------------

import os
import json

import numpy as np
import healpy as hp

import matplotlib
matplotlib.use("Agg")

import distutils.version
mpl_version = distutils.version.LooseVersion(matplotlib.__version__)

from matplotlib import pyplot as plt
plt.rcParams.update({'font.family':'serif', 'text.usetex':True })

from matplotlib.axes import Axes
from matplotlib import text
from matplotlib import ticker
from matplotlib.ticker import Formatter, FixedLocator
from matplotlib import patheffects
from matplotlib.projections.geo import MollweideAxes
from matplotlib.transforms import Transform, Affine2D
from matplotlib.projections import register_projection

### non-standard libraries
from skymap_statistics.detector_cache import detectors ### for detector arms
from skymap_statistics import triangulate

#-------------------------------------------------

axpos = [0.03, 0.03, 0.94, 0.94]

twopi = 2*np.pi
pi2 = 0.5*np.pi

#-------------------------------------------------

### adapted from lalinference (although the source code has since moved...)

if mpl_version >= '1.3.0':
    FixedMollweideAxes = MollweideAxes

elif mpl_version < '1.2.0':
    class FixedMollweideAxes(MollweideAxes):
        """Patched version of matplotlib's Mollweide projection that implements a
        correct inverse transform."""

        name = 'fixed mollweide'

        class FixedMollweideTransform(MollweideAxes.MollweideTransform):

            def inverted(self):
                return FixedMollweideAxes.InvertedFixedMollweideTransform(self._resolution)
            inverted.__doc__ = Transform.inverted.__doc__

        class InvertedFixedMollweideTransform(MollweideAxes.InvertedMollweideTransform):

            def inverted(self):
                return FixedMollweideAxes.FixedMollweideTransform(self._resolution)
            inverted.__doc__ = Transform.inverted.__doc__

            def transform(self, xy):
                x = xy[:, 0:1]
                y = xy[:, 1:2]

                sqrt2 = np.sqrt(2)
                sintheta = y / sqrt2
                with np.errstate(invalid='ignore'):
                    costheta = np.sqrt(1. - 0.5 * y * y)
                longitude = 0.25 * sqrt2 * np.pi * x / costheta
                latitude = np.arcsin(2 / np.pi * (np.arcsin(sintheta) + sintheta * costheta))
                return np.concatenate((longitude, latitude), 1)
            transform.__doc__ = Transform.transform.__doc__

            transform_non_affine = transform

        def _get_core_transform(self, resolution):
            return self.FixedMollweideTransform(resolution)

else:
    class FixedMollweideAxes(MollweideAxes):
        """Patched version of matplotlib's Mollweide projection that implements a
        correct inverse transform."""

        name = 'fixed mollweide'

        class FixedMollweideTransform(MollweideAxes.MollweideTransform):

            def inverted(self):
                return FixedMollweideAxes.InvertedFixedMollweideTransform(self._resolution)
            inverted.__doc__ = Transform.inverted.__doc__

        class InvertedFixedMollweideTransform(MollweideAxes.InvertedMollweideTransform):

            def inverted(self):
                return FixedMollweideAxes.FixedMollweideTransform(self._resolution)
            inverted.__doc__ = Transform.inverted.__doc__

            def transform_non_affine(self, xy):
                x = xy[:, 0:1]
                y = xy[:, 1:2]

                sqrt2 = np.sqrt(2)
                sintheta = y / sqrt2
                with np.errstate(invalid='ignore'):
                    costheta = np.sqrt(1. - 0.5 * y * y)
                longitude = 0.25 * sqrt2 * np.pi * x / costheta
                latitude = np.arcsin(2 / np.pi * (np.arcsin(sintheta) + sintheta * costheta))
                return np.concatenate((longitude, latitude), 1)
            transform_non_affine.__doc__ = Transform.transform_non_affine.__doc__

        def _get_core_transform(self, resolution):
            return self.FixedMollweideTransform(resolution)

class AstroDegreesMollweideAxes(FixedMollweideAxes):
    """Mollweide axes with phi axis flipped and in degrees from 360 to 0
    instead of in degrees from -180 to 180."""

    name = 'astro degrees mollweide'

    def cla(self):
        super(AstroDegreesMollweideAxes, self).cla()
        self.set_xlim(0, 2*np.pi)

    def set_xlim(self, *args, **kwargs):
        Axes.set_xlim(self, 0., 2*np.pi)
        Axes.set_ylim(self, -np.pi / 2.0, np.pi / 2.0)

    def _get_core_transform(self, resolution):
        return Affine2D().translate(-np.pi, 0.) + super(AstroDegreesMollweideAxes, self)._get_core_transform(resolution)

    def set_longitude_grid(self, degrees):
        # Copied from matplotlib.geo.GeoAxes.set_longitude_grid and modified
        super(AstroDegreesMollweideAxes, self).set_longitude_grid(degrees)
        number = (360.0 / degrees) + 1
        self.xaxis.set_major_locator(
            FixedLocator(
                np.linspace(0, 2*np.pi, int(number), True)[1:-1]))

    def _set_lim_and_transforms(self):
        # Copied from matplotlib.geo.GeoAxes._set_lim_and_transforms and modified
        super(AstroDegreesMollweideAxes, self)._set_lim_and_transforms()

        # This is the transform for latitude ticks.
        yaxis_stretch = Affine2D().scale(np.pi * 2.0, 1.0)
        yaxis_space = Affine2D().scale(-1.0, 1.1)
        self._yaxis_transform = \
            yaxis_stretch + \
            self.transData
        yaxis_text_base = \
            yaxis_stretch + \
            self.transProjection + \
            (yaxis_space + \
             self.transAffine + \
             self.transAxes)
        self._yaxis_text1_transform = \
            yaxis_text_base + \
            Affine2D().translate(-8.0, 0.0)
        self._yaxis_text2_transform = \
            yaxis_text_base + \
            Affine2D().translate(8.0, 0.0)

    def _get_affine_transform(self):
        transform = self._get_core_transform(1)
        xscale, _ = transform.transform_point((0, 0))
        _, yscale = transform.transform_point((0, np.pi / 2.0))
        return Affine2D() \
            .scale(0.5 / xscale, 0.5 / yscale) \
            .translate(0.5, 0.5)

class AstroHoursMollweideAxes(AstroDegreesMollweideAxes):
    """Mollweide axes with phi axis flipped and in hours from 24 to 0 instead of
    in degrees from -180 to 180."""

    name = 'astro mollweide'

    class RaFormatter(Formatter):
        # Copied from matplotlib.geo.GeoAxes.ThetaFormatter and modified
        def __init__(self, round_to=1.0):
            self._round_to = round_to

        def __call__(self, x, pos=None):
            hours = (x / np.pi) * 12.
            hours = round(15 * hours / self._round_to) * self._round_to / 15
            return r"%0.0f$^\mathrm{h}$" % hours

    def set_longitude_grid(self, degrees):
        super(AstroHoursMollweideAxes, self).set_longitude_grid(degrees)
        self.xaxis.set_major_formatter(self.RaFormatter(degrees))

### actually register custom projections
register_projection(FixedMollweideAxes)
register_projection(AstroDegreesMollweideAxes)
register_projection(AstroHoursMollweideAxes)

#------------------------

def outline_text(ax):
    """If we are using a new enough version of matplotlib, then
    add a white outline to all text to make it stand out from the background."""
    effects = [patheffects.withStroke(linewidth=2, foreground='w')]
    for artist in ax.findobj(text.Text):
        artist.set_path_effects(effects)

def _healpix_lookup(map, lon, lat, nest=False, dlon=0):
    """Look up the value of a HEALPix map in the pixel containing the point
    with the specified longitude and latitude."""
    nside = hp.npix2nside(len(map))
    return map[hp.ang2pix(nside, 0.5 * np.pi - lat, lon - dlon, nest=nest)]

def healpix_heatmap(post, ax, *args, **kwargs):
    """Produce a heatmap from a HEALPix map."""
    mpl_kwargs = dict(kwargs)
    dlon = mpl_kwargs.pop('dlon', 0)
    nest = mpl_kwargs.pop('nest', False)

    # Set up a regular grid tiling the bounding box of the axes.
    x = np.arange(ax.bbox.x0, ax.bbox.x1 + 0.5, 0.5)
    y = np.arange(ax.bbox.y0, ax.bbox.y1 + 0.5, 0.5)
    xx, yy = np.meshgrid(x, y)

    # Get axis data limits.
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()

    # Retrieve the inverse transform of the current axes (which converts display
    # coodinates to data coordinates).
    itrans = ax.transData.inverted()

    # Get the longitude and latitude of every point in the bounding box.
    lons, lats = itrans.transform(np.column_stack((xx.ravel(), yy.ravel()))).T

    # Create a mask that selects only the pixels that fall inside the map boundary.
    mask = np.isfinite(lons) & np.isfinite(lats) & (lons >= xmin) & (lons <= xmax)
    zz = np.ma.array(np.empty(lons.shape), mask=~mask)

    # Evaluate the function everywhere that the mask is set.
    zz[mask] = _healpix_lookup(post, lons[mask], lats[mask], nest=nest, dlon=dlon)

    # Plot bitmap using imshow.
    aximg = plt.imshow(zz.reshape(xx.shape), aspect=ax.get_aspect(),
        origin='upper', extent=(xmin, xmax, ymax, ymin),
        *args, **kwargs)

    # Hide masked-out values by displaying them in transparent white.
    aximg.cmap.set_bad('w', alpha=0.)

    # Done.
    return aximg

def healpix_contour(post, ax, *args, **kwargs):
    """Produce a contour plot from a HEALPix map."""
    mpl_kwargs = dict(kwargs)
    dlon = mpl_kwargs.pop('dlon', 0)
    nest = mpl_kwargs.pop('nest', False)
    filled = mpl_kwargs.pop('filled', False)

    # Set up a regular grid tiling in right ascension and declination
    x = np.linspace(*ax.get_xlim(), num=500)
    y = np.linspace(*ax.get_ylim(), num=500)
    xx, yy = np.meshgrid(x, y)

    # Evaluate the function everywhere.
    zz = _healpix_lookup(post, xx, yy, nest=nest, dlon=dlon)

    # Add contour plot
    if filled:
        ax = plt.contourf(xx, yy, zz, *args, **kwargs)
    else:
        ax = plt.contour(xx, yy, zz, *args, **kwargs)

    # Done.
    return ax

#-------------------------------------------------

### actual plotting and figure manipulation

def gen_fig_ax( figind, figheight=5, figwidth=9, projection=None, grid=True ):
    '''
    generates figure and axis in the set-up we prefer
    '''
    fig = plt.figure(figind, figsize=(figwidth, figheight) )
    if projection:
        ax = fig.add_axes(axpos, projection=projection)
    else:
        ax = fig.add_axes(axpos)
    ax.grid( grid )

    return fig, ax

def annotate( ax, projection=None, line_of_sight=[], line_of_sight_color='k', line_of_sight_fontsize=8, line_of_sight_markersize=2, zenith=[], zenith_color='k', zenith_markersize=2, zenith_fontsize=8, time_delay=[], time_delay_color='k', time_delay_alpha=1.0, time_delay_linestyle='solid', marker_Dec_RA=[], marker='o', marker_color='k', marker_size=4, marker_edgewidth=1, marker_alpha=1.0, continents=[], continents_color='k', continents_alpha=1.0, constellations=[], constellations_color='k', constellations_alpha=1.0, stars=[], stars_color='k', stars_alpha=1.0, arms=[], arms_color='k', arms_linewidth=1, arms_alpha=1.0, constellation_boundaries=[], constellation_boundaries_color='k', constellation_boundaries_alpha=1.0, constellation_centers=[], constellation_centers_color='k', constellation_centers_alpha=1.0, constellation_centers_fontsize=8 ):
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

        ax.plot( x, y, color=line_of_sight_color, markeredgecolor=line_of_sight_color, marker='o', markersize=line_of_sight_markersize )
        ax.plot( X, Y, color=line_of_sight_color, markeredgecolor=line_of_sight_color, marker='o', markersize=line_of_sight_markersize )

        ax.text( x, y, " %s-%s"%(ifos[1],ifos[0]), ha='left', va='bottom', color=line_of_sight_color, fontsize=line_of_sight_fontsize )
        ax.text( X, Y, " %s-%s"%(ifos[0],ifos[1]), ha='left', va='bottom', color=line_of_sight_color, fontsize=line_of_sight_fontsize )

    ### plot zenith markers
    for ifo, (y,x), (Y,X) in zenith:
        if projection=="mollweide":
            if x > np.pi:
                x -= twopi
            if X > np.pi:
                X -= twopi

        ax.plot( x, y, color=zenith_color, markeredgecolor=zenith_color, marker='s', markersize=zenith_markersize )
        ax.plot( X, Y, color=zenith_color, markeredgecolor=zenith_color, marker='s', markersize=zenith_markersize )

        ax.text( x, y, " "+ifo+"+", ha='left', va='bottom', color=zenith_color, fontsize=zenith_fontsize )
        ax.text( X, Y, " "+ifo+"-", ha='left', va='bottom', color=zenith_color, fontsize=zenith_fontsize )

    ### plot time-delay loci
    for y, x in time_delay:
        ax.plot( x, y, color=time_delay_color, alpha=time_delay_alpha, linestyle=time_delay_linestyle )

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
        ax.plot( verts[:,0], verts[:,1], color=continents_color, linewidth=0.5, alpha=continents_alpha )
        ax.plot( verts[:,0]+twopi, verts[:,1], color=continents_color, linewidth=0.5, alpha=continents_alpha )
        ax.plot( verts[:,0]-twopi, verts[:,1], color=continents_color, linewidth=0.5, alpha=continents_alpha )

    ### add arms
    for x, y in arms:
        ax.plot( x, y, color=arms_color, linewidth=arms_linewidth, alpha=arms_alpha )

    ### add constellations
    for shape in constellations:
        ax.plot( shape[:,0], shape[:,1], color=constellations_color, linewidth=0.5, alpha=constellations_alpha )
        ax.plot( shape[:,0]+twopi, shape[:,1], color=constellations_color, linewidth=0.5, alpha=constellations_alpha )
        ax.plot( shape[:,0]-twopi, shape[:,1], color=constellations_color, linewidth=0.5, alpha=constellations_alpha )

    ### add stars
    for x, y, mag in stars:
        markersize=max(1, 5-mag) ### FIXME: hard coded...bad?
        if projection=="mollweide":
            if x < -np.pi:
                x += twopi
            elif x > np.pi:
                x -= twopi
        elif projection=="astro mollweide":
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
        if projection=="mollweide":
            if x < -np.pi:
                x += twopi
            elif x > np.pi:
                x -= twopi
        elif projection=="astro mollweide":
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

def heatmap( post, ax, color_map='OrRd', colorbar=False, colorbar_label='' ):
    '''
    generate mollweide projection of heatmap with requested annotations
    '''
    healpix_heatmap( post, ax, cmap=plt.get_cmap(color_map) ) ### is this buggy when projection=="mollweide"?
    if colorbar:
        cb = plt.colorbar(orientation='horizontal', fraction=0.15, pad=0.03, shrink=0.8) ### FIXME: hard-coded options are a bit fragile...
        cb.set_label(colorbar_label)

def contour( post, ax, levels=[0.1, 0.5, 0.9], alpha=1.0, colors='b', linewidths=1, filled=False ):
    '''
    generate mollweide projection of contours with requested annotations
    '''
    cpost = np.empty(post.shape)
    indecies = np.argsort(post)[::-1]
    cpost[indecies] = np.cumsum(post[indecies])

    healpix_contour( cpost,
                     ax,
                     levels=levels,
                     alpha=alpha,
                     colors=colors,
                     linewidths=linewidths,
                     filled=filled )

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
            y = pi2 - y
            Y = pi2 - Y
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
            y = pi2 - y
            Y = pi2 - Y
        zenith.append( (ifo, (y,x), (Y,X)) )
    return zenith

def gen_arms( IFOs, coord='C', gps=None, extend=1.0):
    '''
    computes the detector arms and returns them in a friendly format for plotting
    this is assumed to *only* be run for mollweide projections, so I've hard-coded some lengths

    you can extend the length of the arms using "extend"!=1.0
    '''
    arms = []
    for ifo in IFOs:
        ### get these in coord=E
        det = detectors[ifo]

        crnr = det.dr / np.sum(det.dr**2)**0.5
        crnr_y = np.arccos(crnr[2])
        crnr_x = np.arctan2(crnr[1], crnr[0])

        xend = det.dr + det.nx*1e-3 * extend ### extend controls arm length
        xend /= np.sum(xend**2)**0.5
        xend_y = np.arccos(xend[2])
        xend_x = np.arctan2(xend[1], xend[0])

        yend = det.dr + det.ny*1e-3 * extend ### extend controls arm length
        yend /= np.sum(yend**2)**0.5
        yend_y = np.arccos(yend[2])
        yend_x = np.arctan2(yend[1], yend[0])

        ### convert theta->dec
        crnr_y = pi2 - crnr_y 
        xend_y = pi2 - xend_y 
        yend_y = pi2 - yend_y 

        if coord=="C":
            crnr_x = triangulate.rotateRAE2C(crnr_x, gps, noWRAP=True)
            xend_x = triangulate.rotateRAE2C(xend_x, gps, noWRAP=True)
            yend_x = triangulate.rotateRAE2C(yend_x, gps, noWRAP=True)

        arms.append( np.array([(crnr_x, xend_x), (crnr_y, xend_y)]) )
        arms.append( np.array([(crnr_x, yend_x), (crnr_y, yend_y)]) )

    return arms

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
                y = pi2 - y

                x[x>np.pi] -= twopi ### ensure that everything is between -pi and pi

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
                ra -= twopi
            marker_Dec_RA.append( (pi2-dec, ra) )
    else: ### coord=="C"
        for dec, ra in SRCs:
            if degrees:
                dec *= triangulate.deg2rad
                ra *= triangulate.deg2rad
#            if ra > np.pi:
#                ra -= twopi
            marker_Dec_RA.append( (dec, ra) )

    return marker_Dec_RA

def gen_continents( coord='C', gps=None ):
    '''
    extract the little line segments needed to plot continents and ensure they're all in the correct coordinate system
    '''
    geojson_filename = os.path.join(os.path.dirname(__file__), 'ne_simplified_coastline.json')
    file_obj = open(geojson_filename, 'r')
    geojson = json.load(file_obj)
    file_obj.close()

    verts = [np.deg2rad(shape['coordinates']) for shape in geojson['geometries']]

    if coord=='C': ### rotate into Celestial coordinages
        for vert in verts:
            vert[:,0] = triangulate.rotateRAE2C(vert[:,0], gps, noWRAP=True)

    return verts

def gen_constellations( coord='C', gps=None ):
    '''
    extract the constellations from disk and prepare them for plotting
    '''
    json_filename = os.path.join(os.path.dirname(__file__), 'constellationsANDstars.json')
    file_obj = open(json_filename, 'r')
    constjson = json.load(file_obj)
    file_obj.close()

    shapes = []
    for value in constjson['constellations'].values():
        for shape in value:
            shapes.append( np.array([np.array(_) for _ in shape]) )
    shapes = np.array(shapes)

    if coord=='E': ### rotate into E coordinates
        for shape in shapes:
            shape[:,0] = triangulate.rotateRAC2E(shape[:,0], gps, noWRAP=True)

    return shapes

def gen_stars( coord='C', gps=None ):
    '''
    extract the bright stars from disk and prepare them for plotting
    '''
    json_filename = os.path.join(os.path.dirname(__file__), 'constellationsANDstars.json')
    file_obj = open(json_filename, 'r')
    constjson = json.load(file_obj)
    file_obj.close()

    stars = np.array(constjson['stars'])

    if coord=='E': ### rotate into E coordinates
        stars[:,0] = triangulate.rotateRAC2E(stars[:,0], gps)

    return stars

def gen_constellationBoundaries( coord='C', gps=None ):
    '''
    extract constellation boundaries from disk and prepare them for plotting
    '''
    json_filename = os.path.join(os.path.dirname(__file__), 'constellationsANDstars.json')
    file_obj = open(json_filename, 'r')
    constjson = json.load(file_obj)
    file_obj.close()

    boundaries = []
    for edge in constjson['boundaries']:
        boundaries.append( np.array([np.array(edge[0]), np.array(edge[1])]) )
    boundaries = np.array(boundaries)

    if coord=='E':
        boundaries[:,0] = triangulate.rotateRAC2E(boundaries[:,0], gps, noWRAP=True)
    return boundaries

def gen_constellationCenters( coord='C', gps=None ):
    '''
    extract constellatoin centers from disk and prepare them for plotting
    '''
    json_filename = os.path.join(os.path.dirname(__file__), 'constellationsANDstars.json')
    file_obj = open(json_filename, 'r')
    constjson = json.load(file_obj)
    file_obj.close()

    centers = constjson['centers']

    if coord=='E': ### rotate into E coordinates
        for name, (ra, dec) in centers.items():
            centers[name] = (triangulate.rotateRAC2E(ra, gps), dec)

    return [(ra, dec, name) for name, (ra, dec) in centers.items()]
