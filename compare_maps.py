usage = "python compare_maps.py [--options] label1,fits1 label2,fits2 ..."
description = "computes several comparison statistics for the FITs files provided. Computes comparison statistics for all possible pairings (downsampling the finner resolution map if needed)"
author = "R. Essick (reed.essick@ligo.org)"

#==========================================================

import numpy as np
import healpy as hp

import stats

from optparse import OptionParser

#==========================================================
parser = OptionParser(usage=usage, description=description)

parser.add_option("-v", "--verbose", default=False, action="store_true")
parser.add_option("-V", "--Verbose", default=False, action="store_true")

parser.add_option("-d", "--degrees", default=False, action="store_true")

parser.add_option("", "--dMAP", default=False, action="store_true", help="compute the angular separation between maximum a posteriori points")

parser.add_option("", "--fidelity", default=False, action="store_true", help="compute the fidelity between maps")
parser.add_option("", "--joint-entropy", default=False, action="store_true", help="compute joint entropy between maps")

parser.add_option("", "--symKL", default=False, action="store_true", help="compute symmetric KLdivergence between maps")
parser.add_option("", "--mse", default=False, action="store_true", help="compute the mean square error between maps")
parser.add_option("", "--peak-snr", default=False, action="store_true", help="compute peak snr between maps")
parser.add_option("", "--structural-similarity", default=False, action="store_true", help="compute structural similarity between maps")
parser.add_option("", "--pearson", default=False, action="store_true", help="compute pearson correlation coefficient between maps")
parser.add_option("", "--dot", default=False, action="store_true", help="compute normalized dot product between maps")

parser.add_option("-c", "--credible-interval", default=[], type='float', action='append', help='compute the overlap and intersection of the credible intervals reported in the maps')

parser.add_option("", "--graceid", default=[], type="string", action="append", help="will upload annotations to GraceDB events. if used, there must be one graceid per argment. DO NOT USE UNLESS YOU HAVE LALSuite AND PERMISSION TO ANNOTATE GraceDB!")

parser.add_option('', '--gdb-url', default='https://gracedb.ligo.org/api', type='string')
parser.add_option('', '--tag-as-sky-loc', default=False, action='store_true')

parser.add_option("", "--skip-gracedb-upload", default=False, action="store_true")

opts, args = parser.parse_args()

opts.verbose = opts.verbose or opts.Verbose

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

if opts.degrees:
        unit = "deg"
        areaunit = "deg2"
        angle_conversion = 180/np.pi
else:
        unit = "radians"
        areaunit = "stradians"
        angle_conversion = 1.0

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

for ind, label1 in enumerate(labels):
	d1 = maps[label1]
	Post1 = d1['post']
	nside1 = d1['nside']
	if opts.graceid:
		gid1 = d1['graceid']

	for label2 in labels[ind+1:]:

		d2 = maps[label2]
		Post2 = d2['post']
		nside2 = d2['nside']
		if opts.graceid:
			gid2 = d2['graceid']

		print "%s vs %s"%(label1, label2)

		### resample if needed
		nside = min(nside1, nside2)
		pixarea = hp.nside2pixarea(nside, degrees=opts.degrees)

		if nside2 > nside:
			if opts.verbose:
				print "\tresampling %s : %d -> %d"%(label2, nside2, nside)
			post2 = stats.resample(Post2[:], nside)
		else:
			post2 = Post2[:]
		if nside1 > nside:
			if opts.verbose:
				print "\tresampling %s : %d -> %d"%(label1, nside1, nside)
			post1 = stats.resample(Post1[:], nside)
		else:
			post1 = Post1

		pixarea = hp.nside2pixarea( nside, degrees=opts.degrees )

		if opts.verbose:
			print "\tnside : %d"%(nside)
			print "\tpixare : %.6f %s"%(pixarea, areaunit)

		messages = []
	
		### compute statistics
		if opts.dMAP:
			if opts.Verbose:
				print "\t\tdMAP"
			t1, p1 = d1['estang']
			t2, p2 = d2['estang']
			messages.append( "dtheta_MAP : %.5f %s"%(angle_conversion*np.arccos(stats.cos_dtheta(t1, p1, t2, p2, safe=True)), unit) )

		if opts.fidelity:
			if opts.Verbose:
				print "\t\tfidelity"
			messages.append( "fidelity : %.5f"%(stats.fidelity(post1, post2)) )

		if opts.joint_entropy:
			if opts.Verbose:
				print "\t\tjoint_entropy"
			messages.append( "joint entropy : %.5f %s"%(pixarea * 2**(stats.indep_joint_entropy(post1, post2, base=2.0)), areaunit) )

		if opts.symKL:
			if opts.Verbose:
				print "\t\tsymKL"
			messages.append( "symmetric KL divergence : %.5f"%stats.symmetric_KLdivergence(post1, post2) )

		if opts.mse:
			if opts.Verbose:
				print "\t\tmse"
			messages.append( "mean square error : %.5e"%stats.mse(post1, post2) )

		if opts.peak_snr:
			if opts.Verbose:
				print "\t\tpeak_snr"
			messages.append( "peak SNR : (%.5f, %.5f)"%(stats.peak_snr(post1, post2)) )

		if opts.structural_similarity:
			if opts.Verbose:
				print "\t\tstructural_similarity"
			messages.append( "structural similarity : %.5f"%stats.structural_similarity(post1, post2, c1=0.01*np.min(np.mean(post1),np.mean(post2)), c2=0.01*np.min(np.var(post1), np.var(post2)) ) )

		if opts.pearson:
			if opts.Verbose:
				print "\t\tpearson"
			messages.append( "pearson : %.5f"%stats.pearson(post1, post2) )

		if opts.dot:
			if opts.Verbose:
				print "\t\tdot"
			messages.append( "dot : %.5f"%stats.dot(post1, post2) )

		if opts.Verbose:
			print "\t\tCredible Regions"
		for conf, pix1, pix2 in zip(opts.credible_interval, stats.credible_region(post1, opts.credible_interval), stats.credible_region(post2, opts.credible_interval) ):
			if opts.Verbose:
				print "\t\tCR: %.6f"%(conf)
			header = "%.3f %s CR"%(conf*100, "%")

			messages.append( "%s : %s = %.3f %s"%(header, label1, pixarea*len(pix1) , areaunit) )
			messages.append( "%s : %s = %.3f %s"%(header, label2, pixarea*len(pix2) , areaunit) )

			if opts.Verbose:
				print "\t\tgeometric_overlap"
			i, u = stats.geometric_overlap(pix1, pix2, nside=nside, degrees=opts.degrees)
			messages.append( "%s : intersection = %.3f %s"%(header, i, areaunit) )
			messages.append( "%s : union = %.3f %s"%(header, u, areaunit) )

		for message in messages:
			print "\t", message

		if opts.graceid: ### upload to gracedb
			for gid in set([gid1, gid2]): ### if gid's are identical, only report once
				if opts.verbose and (not opts.skip_gracedb_upload):
					print "uploading to GraceID :", gid
				for message in messages:
					if not opts.skip_gracedb_upload:
						if opts.tag_as_sky_loc:
				                        gracedb.writeLog(gid, message="%s,%s : %s"%(label1, label2, message), filename=None, tagname="sky_loc")
						else:
							gracedb.writeLog(gid, message="%s,%s : %s"%(label1, label2, message), filename=None)
			

		

