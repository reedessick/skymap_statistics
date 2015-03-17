usage = "python analyze_maps.py [--options] label1,fits1 label2,fits2 ..."
description = "computes several descriptive statistics about the maps provided"
author = "R. Essick (reed.essick@ligo.org)"

#==========================================================

import numpy as np
import healpy as hp

import stats

from optparse import OptionParser

#==========================================================
parser = OptionParser(usage=usage, description=description)

parser.add_option("-v", "--verbose", default=False, action="store_true")

parser.add_option("-d", "--degrees", default=False, action="store_true")

parser.add_option("", "--pvalue", default=False, type="string", help="\"theta, phi\" in the coordinate system of these maps for which we compute the pvalue (the confidence associated with the minimum credible region that marginally includes this location)")

parser.add_option("-H", "--entropy", default=False, action="store_true", help='computes the entropy of the map')

parser.add_option("-c", "--credible-interval", default=[], type='float', action='append', help='computes the size, max(dtheta), and num/size of disjoint regions for each confidence region. This argument can be supplied multiple times to analyze multiple confidence regions.')
parser.add_option("", "--no-credible-interval-dtheta", default=False, action="store_true", help='does not compute max(dtheta) for each confidence region. This may be desired if computational speed is an issue, because the max(dtheta) algorithm scales as Npix^2.')
parser.add_option("", "--no-disjoint-regions", default=False, action="store_true", help="does not compute number,size of disjoint regions for each confidence region.")

parser.add_option("", "--graceid", default=[], type="string", action="append", help="will upload annotations to GraceDB events. if used, there must be one graceid per argument. DO NOT USE UNLESS YOU HAVE LALSuite AND PERMISSION TO ANNOTATE GraceDB!")

parser.add_option('', '--gdb-url', default='https://gracedb.ligo.org/api', type='string')
parser.add_option('', '--tag-as-sky-loc', default=False, action='store_true')

opts, args = parser.parse_args()

if opts.graceid:
	from ligo.gracedb.rest import GraceDb
	gracedb = gracedb = GraceDb(opts.gdb_url)

if opts.graceid and len(opts.gracedb)!=len(args):
        raise ValueError("when supplying --graceid, you must supply the same number of graceid entries and fits files")

if opts.pvalue:
	theta, phi = [float(l) for l in opts.pvalue.split(",")]

opts.credible_interval = sorted(opts.credible_interval)

#==========================================================
if opts.degrees:
	unit = "deg"
	areaunit = "deg2"
	angle_conversion = 180/np.pi
	theta /= angle_conversion
	phi /= angle_conversion
else:
	unit = "radians"
	areaunit = "stradians"
	angle_conversion = 1.0

#==========================================================

for ind, arg in enumerate(args):
	label, fits = arg.split(',')
	messages = []

	print label

	if opts.verbose:
		print "\treading map from %s"%(fits)
	post, header = hp.read_map( fits, h=True )
	NEST = (dict(header)['ORDERING']=='NEST')
	npix = len(post)
	nside = hp.npix2nside(npix)

	print "\tnside=%d"%(nside)

	pixarea = hp.nside2pixarea(nside, degrees=opts.degrees)

	### compute statistics and report them
	if opts.pvalue:
		pvalue = stats.p_value(post, theta, phi, nside=nside)
		messages.append( "cdf(%s) = %.3f %s"%(opts.pvalue, pvalue*100, "%") )
		
	# entropy -> size
	if opts.entropy:
		entropy = pixarea * 2**(stats.entropy(post, base=2.0))
		messages.append( "entropy = %.3f %s"%(entropy, areaunit) )

	# CR -> size, max(dtheta)
	cr = {}
	for CR, conf in zip(stats.credible_region(post, opts.credible_interval), opts.credible_interval):
		header = "%.3f %s CR"%(conf*100, "%")
		size = pixarea*len(CR)
		messages.append( "%s: size= %.3f %s"%(header, size, areaunit) )

		if not opts.no_credible_interval_dtheta:
			max_dtheta = angle_conversion*np.arccos(stats.min_all_cos_dtheta(CR, nside, nest=NEST))
			messages.append( "%s: max(dtheta) = %.3f %s"%(header, max_dtheta, unit) )

		if not opts.no_disjoint_regions:
			sizes = sorted([len(_)*pixarea for _ in stats.__into_modes(nside, CR, nest=NEST)])
			messages.append( "%s: disjoint regions : (%s) %s"%(header, ", ".join(["%.3f"%x for x in sizes]), areaunit ) )

	for message in messages:
		print "\t", message

	if opts.graceid: ### upload to GraceDB
		gid = opts.graceid[ind] 
		for message in messages:
			if opts.tag_as_sky_loc:
				gracedb.writeLog(gid, message="%s : %s"%(label, message), filename=None, tagname="sky_loc")
			else:
				gracedb.writeLog(gid, message="%s : %s"%(label, message), filename=None)

		


