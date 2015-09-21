usage = "python compare_maps.py [--options] label1,fits1 label2,fits2 ..."
description = "overlays skymaps on a figure for visualizatoin"
author = "R. Essick (reed.essick@ligo.org)"

#==========================================================

try:
	from lalinference import plot as lalinf_plot
except:
	raise StandardError("Could not import lalinference.plot")

import numpy as np
import healpy as hp

import stats

from mpl_toolkits.basemap import Basemap
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import matplotlib.cm as cm

from optparse import OptionParser

#==========================================================
parser = OptionParser(usage=usage, description=description)

parser.add_option("-v", "--verbose", default=False, action="store_true")

parser.add_option("-c", "--credible-interval", default=[], type='float', action='append', help='compute the overlap and intersection of the credible intervals reported in the maps')

parser.add_option("", "--graceid", default=[], type="string", action="append", help="will upload annotations to GraceDB events. if used, there must be one graceid per argment. DO NOT USE UNLESS YOU HAVE LALSuite AND PERMISSION TO ANNOTATE GraceDB!")

parser.add_option('', '--gdb-url', default='https://gracedb.ligo.org/api', type='string')
parser.add_option('', '--tag-as-sky-loc', default=False, action='store_true')

parser.add_option("", "--skip-gracedb-upload", default=False, action="store_true")

opts, args = parser.parse_args()

if opts.graceid:
        from ligo.gracedb.rest import GraceDb
        gracedb = gracedb = GraceDb(opts.gdb_url)

if opts.graceid and len(opts.graceid)!=len(args):
        raise ValueError("when supplying --graceid, you must supply the same number of graceid entries and fits files")

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

        d['post'] = post
        d['npix'] = npix
        d['nside'] = nside
	d['estang'] = stats.estang(post, nside=nside)

#=================================================
### iterate through pairs and compute statistics

figind = 0
for ind, label1 in enumerate(labels):
	d1 = maps[label1]
	post1 = d1['post']
	nside1 = d1['nside']
	if opts.graceid:
		gid1 = d1['graceid']

	for label2 in labels[ind+1:]:

		d2 = maps[label2]
		post2 = d2['post']
		nside2 = d2['nside']
		if opts.graceid:
			gid2 = d2['graceid']

		print "%s vs %s"%(label1, label2)
		
		fig = plt.figure(figind)

		lalinf_plot.contour( post1 )
#		lalinf_plot.contour( post2 )	

		figname = "%s-%s.png"%(label1, label2)
		if opts.verbose:
			print "\t", figname
		fig.savefig( figname)

		figind += 1


