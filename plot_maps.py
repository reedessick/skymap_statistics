#!/usr/bin/python
usage = "python compare_maps.py [--options] label1,fits1 label2,fits2 ..."
description = "plot skymaps on a figure for visualization. Basically a wrapper for lalinference.plot.healpix_heatmap"
author = "R. Essick (reed.essick@ligo.org)"

#==========================================================

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
try:
	from lalinference import plot as lalinf_plot
	from lalinference import cmap
except:
	raise StandardError("Could not import lalinference.plot")

import numpy as np
import healpy as hp

import stats
import triangulate

from optparse import OptionParser

#==========================================================

colors = ['b', 'r', 'g', 'm', 'c', 'k', 'y']

#==========================================================
parser = OptionParser(usage=usage, description=description)

parser.add_option("-v", "--verbose", default=False, action="store_true")

parser.add_option("", "--stack-posteriors", default=False, action="store_true")

parser.add_option("-l", "--logarithmic", default=False, action="store_true")

parser.add_option("-H", "--figheight", default=5, type="float")
parser.add_option("-W", "--figwidth", default=9, type="float")

parser.add_option("-o", "--output-dir", default=".", type="string")

parser.add_option("", "--graceid", default=[], type="string", action="append", help="will upload annotations to GraceDB events. if used, there must be one graceid per argment. DO NOT USE UNLESS YOU HAVE LALSuite AND PERMISSION TO ANNOTATE GraceDB! Not yet implemented.")

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
parser.add_option("", "--zenith", default=[], action="append", type="string", help="eg: H")
parser.add_option("", "--gps", default=None, type="float", help="must be specified if --line-of-sight or --zenith is used")
parser.add_option("", "--coord", default="C", type="string", help="coordinate system of the maps. Default is celestial (C), but we also know Earth-Fixed (E)")

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

if (opts.line_of_sight or opts.zenith) and (opts.coord!="E") and (opts.gps==None):
    opts.gps = float(raw_input("gps = "))

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

#==========================================================
### load posteriors from fits files

figind = 0
if opts.stack_posteriors:
	stack_fig = plt.figure(figind, figsize=(opts.figwidth, opts.figheight) )
	if opts.projection:
        	stack_ax = plt.subplot(111, projection=opts.projection)
	else:
        	stack_ax = plt.subplot(111)
	stack_ax.grid( True )
	figind += 1

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

	fig = plt.figure( figind, figsize=(opts.figwidth, opts.figheight) )
	if opts.projection:
		ax = plt.subplot(111, projection=opts.projection)
	else:
		ax = plt.subplot(111)
	ax.grid( True )

	lalinf_plot.healpix_heatmap( post, cmap=plt.get_cmap(opts.color_map) )

	for ifos, (y,x), (Y,X) in line_of_sight:
		if x > np.pi:
	        	ax.plot( x-2*np.pi, y, color='k', marker='o', markersize=2 )
        	        ax.text( x-2*np.pi, y, " %s-%s"%(ifos[1],ifos[0]), ha='left', va='bottom' )
		else:
	        	ax.plot( x, y, color='k', marker='o', markersize=2 )
        	        ax.text( x, y, " %s-%s"%(ifos[1],ifos[0]), ha='left', va='bottom' )
		if X > np.pi:	
        	        ax.plot( X-2*np.pi, Y, color='k', marker='o', markersize=2 )
	                ax.text( X-2*np.pi, Y, " %s-%s"%(ifos[0],ifos[1]), ha='left', va='bottom' )
		else:
                	ax.plot( X, Y, color='k', marker='o', markersize=2 )
	                ax.text( X, Y, " %s-%s"%(ifos[0],ifos[1]), ha='left', va='bottom' )

	for ifo, (y,x), (Y,X) in zenith:
		if x > np.pi:
	        	ax.plot( x-2*np.pi, y, color='k', marker='s', markersize=2 )
                	ax.text( x-2*np.pi, y, " "+ifo+"+", ha='left', va='bottom' )
		else:
	        	ax.plot( x, y, color='k', marker='s', markersize=2 )
        	        ax.text( x, y, " "+ifo+"+", ha='left', va='bottom' )
		if X > np.pi:
	        	ax.plot( X-2*np.pi, Y, color='k', marker='s', markersize=2 )
                	ax.text( X-2*np.pi, Y, " "+ifo+"-", ha='left', va='bottom' )
		else:
	        	ax.plot( X, Y, color='k', marker='s', markersize=2 )
        	        ax.text( X, Y, " "+ifo+"-", ha='left', va='bottom' )

	if opts.transparent:
		fig.patch.set_alpha(0.)
		ax.patch.set_alpha(0.)
		ax.set_alpha(0.)

	for figtype in opts.figtype:
	    	figname = "%s/%s%s.%s"%(opts.output_dir, label, opts.tag, figtype)
		if opts.verbose:
			print "\t", figname
		plt.savefig( figname, dpi=opts.dpi )
	plt.close( fig )

        if opts.stack_posteriors:
		plt.sca( stack_ax )
#		lalinf_plot.healpix_heatmap( post, cmap=plt.get_cmap(opts.color_map) )
		lalinf_plot.healpix_contour( cpost, levels=[0.1, 0.5, 0.9], alpha=0.75, label=label, colors=colors[(figind-1)%len(colors)] )
		stack_fig.text(0.01, 0.99-0.05*(figind-1), label, color=colors[(figind-1)%len(colors)], ha='left', va='top')

	figind += 1

if opts.stack_posteriors:
	plt.figure( 0 )
	plt.sca( stack_ax )

        for ifos, (y,x), (Y,X) in line_of_sight:
                if x > np.pi:
                        stack_ax.plot( x-2*np.pi, y, color='k', marker='o', markersize=2 )
                        stack_ax.text( x-2*np.pi, y, " %s-%s"%(ifos[1],ifos[0]), ha='left', va='bottom' )
                else:
                        stack_ax.plot( x, y, color='k', marker='o', markersize=2 )
                        stack_ax.text( x, y, " %s-%s"%(ifos[1],ifos[0]), ha='left', va='bottom' )
                if X > np.pi:
                        stack_ax.plot( X-2*np.pi, Y, color='k', marker='o', markersize=2 )
                        stack_ax.text( X-2*np.pi, Y, " %s-%s"%(ifos[0],ifos[1]), ha='left', va='bottom' )
                else:
                        stack_ax.plot( X, Y, color='k', marker='o', markersize=2 )
                        stack_ax.text( X, Y, " %s-%s"%(ifos[0],ifos[1]), ha='left', va='bottom' )

        for ifo, (y,x), (Y,X) in zenith:
                if x > np.pi:
                        stack_ax.plot( x-2*np.pi, y, color='k', marker='s', markersize=2 )
                        stack_ax.text( x-2*np.pi, y, " "+ifo+"+", ha='left', va='bottom' )
                else:
                        stack_ax.plot( x, y, color='k', marker='s', markersize=2 )
                        stack_ax.text( x, y, " "+ifo+"+", ha='left', va='bottom' )
                if X > np.pi:
                        stack_ax.plot( X-2*np.pi, Y, color='k', marker='s', markersize=2 )
                        stack_ax.text( X-2*np.pi, Y, " "+ifo+"-", ha='left', va='bottom' )
                else:
                        stack_ax.plot( X, Y, color='k', marker='s', markersize=2 )
                        stack_ax.text( X, Y, " "+ifo+"-", ha='left', va='bottom' )

        if opts.transparent:
                stack_fig.patch.set_alpha(0.)
                stack_ax.patch.set_alpha(0.)
                stack_ax.set_alpha(0.)

        for figtype in opts.figtype:
                figname = "%s/stackedPosterior%s.%s"%(opts.output_dir, opts.tag, figtype)
                if opts.verbose:
                        print "\t", figname
                plt.savefig( figname, dpi=opts.dpi )
        plt.close( stack_fig )

