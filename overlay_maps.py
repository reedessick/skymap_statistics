#!/usr/bin/python
usage = "python overlay_maps.py [--options] label1,fits1 label2,fits2 ..."
description = "overlays skymaps on a figure for visualization. Basically a wrapper for lalinference.plot.healpix_contour"
author = "R. Essick (reed.essick@ligo.org)"

#==========================================================

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
plt.rcParams.update( {"text.usetex":True} )
try:
	from lalinference import plot as lalinf_plot
except:
	raise StandardError("Could not import lalinference.plot")

import numpy as np
import healpy as hp

import stats
import triangulate

from optparse import OptionParser

#==========================================================
parser = OptionParser(usage=usage, description=description)

parser.add_option("-v", "--verbose", default=False, action="store_true")

parser.add_option("", "--background", default=None, type="string", help="a FITS file to plot in the background of the stacked plot")

parser.add_option("-c", "--credible-interval", default=[], type='float', action='append', help='compute the overlap and intersection of the credible intervals reported in the maps')

parser.add_option("-H", "--figheight", default=5, type="float")
parser.add_option("-W", "--figwidth", default=9, type="float")

parser.add_option("-o", "--output-dir", default=".", type="string")

parser.add_option("", "--graceid", default=[], type="string", action="append", help="will upload annotations to GraceDB events. if used, there must be one graceid per argment. DO NOT USE UNLESS YOU HAVE LALSuite AND PERMISSION TO ANNOTATE GraceDB! Not yet implemented")

parser.add_option('', '--gdb-url', default='https://gracedb.ligo.org/api', type='string')
parser.add_option('', '--tag-as-sky-loc', default=False, action='store_true')

parser.add_option("", "--skip-gracedb-upload", default=False, action="store_true")

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

parser.add_option("", "--gps", default=None, type="float", help="must be specified if --line-of-sight or --zenith is used")
parser.add_option("", "--coord", default="C", type="string", help="coordinate system of the maps. Default is celestial (C), but we also know Earth-Fixed (E)")

parser.add_option("", "--outline-labels", default=False, action="store_true", help="put a white outline around axis labels")

opts, args = parser.parse_args()

if not opts.figtype:
    opts.figtype.append( "png" )

if opts.tag:
	opts.tag = "_%s"%opts.tag

if opts.graceid:
        from ligo.gracedb.rest import GraceDb
        gracedb = gracedb = GraceDb(opts.gdb_url)

if opts.graceid and len(opts.graceid)!=len(args):
        raise ValueError("when supplying --graceid, you must supply the same number of graceid entries and fits files")

if not opts.credible_interval:
	opts.credible_interval = [0.50, 0.90]
else:
	opts.credible_interval = sorted(opts.credible_interval)

maps = {}
if opts.graceid:
	for gid, arg in zip(opts.graceid, args):
		label, fits = arg.split(",")
		maps[label] = {"fits":fits, "graceid":gid}
else:
	for arg in args:
		label, fits = arg.split(",")
		maps[label] = {"fits":fits}

labels = sorted(maps.keys())

if (opts.line_of_sight or opts.zenith or (opts.time_delay and opts.time_delay_Dec_RA)) and (opts.gps==None):
    opts.gps = float(raw_input("gps = "))

#==========================================================
### load posteriors from fits files

for label in labels:
        d = maps[label]
        fits = d['fits']
        if opts.verbose:
                print "reading map from", fits
	post, header = hp.read_map( fits, h=True )
        npix = len(post)
        if (dict(header)['ORDERING']=='NEST'): ### convert to RING ordering
		post = hp.nest2ring(nside, post)
        nside = hp.npix2nside(npix)
        if opts.verbose:
                print "\tnside=%d"%nside

	cpost = np.empty(post.shape)
	indecies = np.argsort(post)[::-1]
	cpost[indecies] = np.cumsum(post[indecies])

        d['post'] = post
	d['cpost'] = cpost
        d['npix'] = npix
        d['nside'] = nside
	d['estang'] = stats.estang(post, nside=nside)

#=================================================

### figure out positions for line-of-sight and zenith markers
if opts.line_of_sight:
    line_of_sight = []
    for ifos in opts.line_of_sight:
        y, x = triangulate.line_of_sight(ifos[1], ifos[0], coord=opts.coord, tgeocent=opts.gps, degrees=False)
        X, Y = triangulate.antipode( x, y, coord=opts.coord, degrees=False)
        if opts.coord=="E": ### convert theta->dec
            y = 0.5*np.pi - y
            Y = 0.5*np.pi - Y
        line_of_sight.append( (ifos, (y,x), (Y,X)) )
else:
    line_of_sight = []

if opts.zenith:
    zenith = []
    for ifo in opts.zenith:
        y, x = triangulate.overhead(ifo, coord=opts.coord, tgeocent=opts.gps, degrees=False)
        X, Y = triangulate.antipode( x, y, coord=opts.coord, degrees=False)
        if opts.coord=="E": ### convert theta->dec
            y = 0.5*np.pi - y
            Y = 0.5*np.pi - Y
        zenith.append( (ifo, (y,x), (Y,X)) )
else:
    zenith = []

### figure out points for time_delay
if opts.time_delay:
    time_delay = []
    for dec, ra in opts.time_delay_Dec_RA:
        if opts.time_delay_degrees:
            dec *= triangulate.deg2rad
            ra *= triangulate.deg2rad
        for ifos in opts.time_delay:
            dt = triangulate.time_delay( dec, ra, ifos[1], ifos[0], coord=opts.coord, tgeocent=opts.gps, degrees=False )
            y, x = triangulate.time_delay_locus( dt, ifos[1], ifos[0], coord=opts.coord, tgeocent=opts.gps, degrees=False )

            if opts.coord=="E": ### convert theta-> dec
                y = 0.5*np.pi - y

            x[x>np.pi] -= 2*np.pi ### ensure that everything is between -pi and pi

            ### find big jumps in azimuthal angle and split up the plotting jobs
            d = np.concatenate( ([0],np.nonzero(np.abs(x[1:]-x[:-1])>np.pi)[0]+1,[len(x)]) )
            for istart, iend in zip(d[:-1], d[1:]):
                time_delay.append( (y[istart:iend], x[istart:iend]) )

else:
    time_delay = []

### figure out points for markers
if opts.marker_Dec_RA:
    if opts.coord=="E":
        marker_Dec_RA=[]
        for dec, ra in opts.marker_Dec_RA:
            if opts.marker_degrees:
                dec *= triangulate.deg2rad
                ra *= triangulate.deg2rad
            if ra > np.pi:
                ra -= 2*np.pi
            marker_Dec_RA.append( (0.5*np.pi-dec, ra) )
    else: ### opts.coord=="C"
        marker_Dec_RA = []
        for dec, ra in opts.marker_Dec_RA:
            if opts.marker_degrees:
                dec *= triangulate.deg2rad
                ra *= triangulate.deg2rad
            if ra > np.pi:
                ra -= 2*np.pi
            marker_Dec_RA.append( (dec, ra) )
else:
    marker_Dec_RA = []

#=================================================

### iterate through and plot

figind = 0
if opts.background:
    bkgnd = hp.read_map(opts.background)
for ind, label1 in enumerate(labels):
	d1 = maps[label1]
	post1 = d1['post']
	cpost1 = d1['cpost']
	nside1 = d1['nside']
	if opts.graceid:
		gid1 = d1['graceid']

	for label2 in labels[ind+1:]:

		d2 = maps[label2]
		post2 = d2['post']
		cpost2 = d2['cpost']
		nside2 = d2['nside']
		if opts.graceid:
			gid2 = d2['graceid']

		print "%s vs %s"%(label1, label2)
		
		fig = plt.figure( figind, figsize=(opts.figwidth, opts.figheight) )
		if opts.projection:
			ax = plt.subplot(111, projection=opts.projection)
		else:
			ax = plt.subplot(111)
		ax.grid( True )

                if opts.background:
                    lalinf_plot.healpix_heatmap( bkgnd, cmap=plt.get_cmap(opts.color_map) )
		c1 = lalinf_plot.healpix_contour( cpost1, levels=opts.credible_interval, colors='b', alpha=0.75, label=label1 )
		c2 = lalinf_plot.healpix_contour( cpost2, levels=opts.credible_interval, colors='r', alpha=0.75, label=label2 )
#		for c in c2.collections:
#			c.set_linestyle('dashed')

#		ax.legend(loc='upper right')
		fig.text(0.1, 0.9, label1.replace("_","\_"), color='b', ha='center', va='center')
		fig.text(0.9, 0.9, label2.replace("_","\_"), color='r', ha='center', va='center')

                for ifos, (y,x), (Y,X) in line_of_sight:
                    if x > np.pi:
                        ax.plot( x-2*np.pi, y, color=opts.line_of_sight_color, markeredgecolor=opts.line_of_sight_color, marker='o', markersize=2 )
                        ax.text( x-2*np.pi, y, " %s-%s"%(ifos[1],ifos[0]), ha='left', va='bottom', color=opts.line_of_sight_color )
                    else:
                        ax.plot( x, y, color=opts.line_of_sight_color, markeredgecolor=opts.line_of_sight_color, marker='o', markersize=2 )
                        ax.text( x, y, " %s-%s"%(ifos[1],ifos[0]), ha='left', va='bottom', color=opts.line_of_sight_color )
                    if X > np.pi:
                        ax.plot( X-2*np.pi, Y, color=opts.line_of_sight_color, markeredgecolor=opts.line_of_sight_color, marker='o', markersize=2 )
                        ax.text( X-2*np.pi, Y, " %s-%s"%(ifos[0],ifos[1]), ha='left', va='bottom', color=opts.line_of_sight_color )
                    else:
                        ax.plot( X, Y, color=opts.line_of_sight_color, markeredgecolor=opts.line_of_sight_color, marker='o', markersize=2 )
                        ax.text( X, Y, " %s-%s"%(ifos[0],ifos[1]), ha='left', va='bottom', color=opts.line_of_sight_color )

                for ifo, (y,x), (Y,X) in zenith:
                    if x > np.pi:
                        ax.plot( x-2*np.pi, y, color=opts.zenith_color, markeredgecolor=opts.zenith_color, marker='s', markersize=2 )
                        ax.text( x-2*np.pi, y, " "+ifo+"+", ha='left', va='bottom', color=opts.zenith_color )
                    else:
                        ax.plot( x, y, color=opts.zenith_color, markeredgecolor=opts.zenith_color, marker='s', markersize=2 )
                        ax.text( x, y, " "+ifo+"+", ha='left', va='bottom', color=opts.zenith_color )
                    if X > np.pi:
                        ax.plot( X-2*np.pi, Y, color=opts.zenith_color, markeredgecolor=opts.zenith_color, marker='s', markersize=2 )
                        ax.text( X-2*np.pi, Y, " "+ifo+"-", ha='left', va='bottom', color=opts.zenith_color )
                    else:
                        ax.plot( X, Y, color=opts.zenith_color, markeredgecolor=opts.zenith_color, marker='s', markersize=2 )
                        ax.text( X, Y, " "+ifo+"-", ha='left', va='bottom', color=opts.zenith_color )

                for y, x in time_delay:
                    ax.plot( x, y, color=opts.time_delay_color, alpha=opts.time_delay_alpha )

                for dec, ra in marker_Dec_RA:
                    ax.plot( ra, dec, linestyle='none', marker='o', markerfacecolor='none', markeredgecolor=opts.marker_color, markersize=4, alpha=opts.marker_alpha )

                if opts.transparent:
                    fig.patch.set_alpha(0.)
                    ax.patch.set_alpha(0.)
                    ax.set_alpha(0.)

                if opts.outline_labels:
		    lalinf_plot.outline_text(stack_ax)


		for figtype in opts.figtype:
  			figname = "%s/%s-%s%s.%s"%(opts.output_dir, label1, label2, opts.tag, figtype)
			if opts.verbose:
				print "\t", figname
			plt.savefig( figname, dpi=opts.dpi )
			plt.close( fig )

		figind += 1
