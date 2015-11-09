usage = "python compare_maps.py [--options] label1,fits1 label2,fits2 ..."
description = "plot skymaps on a figure for visualization. Basically a wrapper for lalinference.plot.healpix_heatmap"
author = "R. Essick (reed.essick@ligo.org)"

#==========================================================

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
try:
	from lalinference import plot as lalinf_plot
except:
	raise StandardError("Could not import lalinference.plot")

import numpy as np
import healpy as hp

import stats

from optparse import OptionParser

#==========================================================
parser = OptionParser(usage=usage, description=description)

parser.add_option("-v", "--verbose", default=False, action="store_true")

parser.add_option("-l", "--logarithmic", default=False, action="store_true")

parser.add_option("-H", "--figheight", default=5, type="float")
parser.add_option("-W", "--figwidth", default=9, type="float")

parser.add_option("-o", "--output-dir", default=".", type="string")

parser.add_option("", "--graceid", default=[], type="string", action="append", help="will upload annotations to GraceDB events. if used, there must be one graceid per argment. DO NOT USE UNLESS YOU HAVE LALSuite AND PERMISSION TO ANNOTATE GraceDB!")

parser.add_option('', '--gdb-url', default='https://gracedb.ligo.org/api', type='string')
parser.add_option('', '--tag-as-sky-loc', default=False, action='store_true')

parser.add_option("", "--skip-gracedb-upload", default=False, action="store_true")

parser.add_option("-p", "--projection", default="astro mollweide", type="string", help="either \"mollweide\", \"astro mollweide\"")

parser.add_option("-t", "--tag", default="", type="string")

opts, args = parser.parse_args()

if opts.tag:
	opts.tag = "_%s"%opts.tag

if opts.graceid:
        from ligo.gracedb.rest import GraceDb
        gracedb = gracedb = GraceDb(opts.gdb_url)

if opts.graceid and len(opts.graceid)!=len(args):
        raise ValueError("when supplying --graceid, you must supply the same number of graceid entries and fits files")

#if not opts.credible_interval:
#	opts.credible_interval = [0.50, 0.90]
#else:
#	opts.credible_interval = sorted(opts.credible_interval)

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

#==========================================================

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

	raise StandardError("WRITE PLOTTING ROUTINES")



"""
#=================================================
### iterate through pairs and compute statistics

figind = 0
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

		c1 = lalinf_plot.healpix_contour( cpost1, levels=opts.credible_interval, colors='b', alpha=0.75, label=label1 )
		c2 = lalinf_plot.healpix_contour( cpost2, levels=opts.credible_interval, colors='r', alpha=0.75, label=label2 )
#		for c in c2.collections:
#			c.set_linestyle('dashed')

#		ax.legend(loc='upper right')
		fig.text(0.1, 0.9, label1, color='b', ha='center', va='center')
		fig.text(0.9, 0.9, label2, color='r', ha='center', va='center')

		figname = "%s/%s-%s%s.png"%(opts.output_dir, label1, label2, opts.tag)
		if opts.verbose:
			print "\t", figname
		plt.savefig( figname)
#		fig.savefig( figname)

		figind += 1
"""

