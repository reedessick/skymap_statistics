usage = "grinch_executable.py [--options] graceid"
description = "manages the interface between GraceDB and the command line inputs for skymap quantification code"
author = "R. Essick (reed.essick@ligo.org), A. Urban (aurban@uwm.edu)"

#=================================================

import os

from math import ceil, floor

from collections import defaultdict
import json

import subprocess

from ligo.gracedb.rest import GraceDb

import ConfigParser
from optparse import OptionParser

#=================================================

parser = OptionParser(usage=usage, description=description)

parser.add_option("-v", "--verbose", default=False, action="store_true")

parser.add_option("-c", "--config", default="grinch_config.ini", type="string")

parser.add_option("", "--analyze", default=False, action="store_true")
parser.add_option("", "--compare", default=False, action="store_true")

parser.add_option("", "--skip-gracedb-upload", default=False, action="store_true")
parser.add_option("", "--keep-maps", default=False, action="store_true")

opts, args = parser.parse_args()

if len(args)!= 1:
	raise ValueError("please supply exactly one GraceID as an argument!")
graceid = args[0]

if not (opts.analyze or opts.compare):
	if opts.verbose:
		print "nothing to do... please supply at least one of \"--analyze\", \"--compare\""
	import sys
	sys.exit(0)

#=================================================

if opts.verbose: 
	print "reading config from", opts.config

config = ConfigParser.SafeConfigParser()
config.read(opts.config)

### used to upload summary files to GraceDB
gdburl = config.get('general', 'gdb_url')

### instantiate gracedb interface
if opts.verbose: 
	print "instantiating connection to GraceDb :", gdburl
gracedb = GraceDb( gdburl )

#=================================================

if opts.verbose:
	print "setting up command lines"

### get universal options for all executables
universal_cmd = " ".join(" -c %s "%(conf) for conf in config.get('general', 'credible_interval').split() )
if config.getboolean('general', 'degrees'):
	universal_cmd += " --degrees "
universal_cmd += " --gdb-url %s "%gdburl

### build analyze command line
analyze_cmd = config.get('executables','analyze_maps') + universal_cmd
if opts.verbose:
	analyze_cmd += " --verbose "
if config.has_option('analyze_maps', 'pvalue'):
	analyze_cmd += " --pvalue %s"%(config.get('analyze_maps','pvalue'))
if config.getboolean('analyze_maps', 'entropy'):
	analyze_cmd += " --entropy "
if config.getboolean('analyze_maps', 'no_credible_interval_dtheta'):
	analyze_cmd += " --no-credible-interval-dtheta "
if config.getboolean('analyze_maps', 'no_credible_interval_disjoint_regions'):
	analyze_cmd += " --no-credible-interval-disjoint-regions "
if config.getboolean('analyze_maps', 'tag_as_sky_loc'):
	analyze_cmd += " --tag-as-sky-loc "

### build compare command line 
compare_cmd = config.get('executables','compare_maps') + universal_cmd
if opts.verbose:
	compare_cmd += " --verbose "
if config.getboolean('compare_maps', 'dMAP'):
	compare_cmd += " --dMAP "
if config.getboolean('compare_maps', 'fidelity'):
	compare_cmd += " --fidelity "
if config.getboolean('compare_maps', 'joint_entropy'):
	compare_cmd += " --joint-entropy "
if config.getboolean('compare_maps', 'symKL'):
	compare_cmd += " --symKL "
if config.getboolean('compare_maps', 'mse'):
	compare_cmd += " --mse "
if config.getboolean('compare_maps', 'peak_snr'):
	compare_cmd += " --peak-snr "
if config.getboolean('compare_maps', 'structural_similarity'):
	compare_cmd += " --structural-similarity "
if config.getboolean('compare_maps', 'pearson'):
	compare_cmd += " --pearson "
if config.getboolean('compare_maps', 'dot'):
	compare_cmd += " --dot "
if config.getboolean('compare_maps', 'tag_as_sky_loc'):
	compare_cmd += " --tag-as-sky-loc "


if opts.skip_gracedb_upload:
	analyze_cmd += " --skip-gracedb-upload"
	compare_cmd += " --skip-gracedb-upload"


#=================================================

names = [] ### list of files that need to be cleaned up after we're done

gdb_entries = { graceid : gracedb.event( graceid ).json() }

#=================================================
### analyze_maps
### only run on the graceid supplied through arguments!
#=================================================
if opts.analyze:
	if opts.verbose:
		print "determining which fits files for %s need to be analyzed"%graceid

	fits = [ _ for _ in gracedb.files( graceid ).json().keys() if _.endswith("fits") or _.endswith("fits.gz") ]
	analyzed = set()

	### get log entries and sort them by relevant fits files
	for log in gracedb.logs( graceid ).json()['log']:
		for fit in fits:
			if "%s-%s : "%(graceid, fit) in log['comment']:
				analyzed.add( fit )

	### download files to local directory
	for fit in [ fit for fit in fits if fit not in analyzed ]:
		name = "%s-%s"%(graceid, fit)
		if opts.verbose:
			print "\tdownloading %s from %s to %s"%(fit, graceid, name)
		file_obj = open(name, "w")
		file_obj.write( gracedb.files( graceid, fit ).read() )
		file_obj.close()
		names.append( name )

		### add them to analyze command
		analyze_cmd += " %s,%s --graceid %s "%(name, name, graceid)
		
	### launch analyze command!
	if opts.verbose:
		print analyze_cmd
	proc = subprocess.Popen(analyze_cmd.split())
	proc.wait() # block!
	if opts.verbose:
		print "\tdone!"

#=================================================
### compare_maps
### run over neighbors!
#=================================================	

if opts.compare:
	### find neighbors
	t = float(gdb_entries[graceid]['gpstime'])
	p_dt = config.getfloat('general', '+dt')
	m_dt = config.getfloat('general', '-dt')

	if opts.verbose:
		print "finding neighbors within [%.5f, %.5f]"%(t+m_dt, t+p_dt)

	gdb_entries.update( dict( (gdb_entry['graceid'], gdb_entry) for gdb_entry in gracedb.events( "%d .. %d"%(floor(t+m_dt), ceil(t+p_dt)) ) if gdb_entry['graceid']!=graceid ) )
	graceids = sorted(gdb_entries.keys()) ### sorting assures the correct ordering within dictionary keys

	if opts.verbose:
		print "\tfound %d neighbors"%(len(graceids)-1)
		for gid in graceids:
			if gid != graceid:
				print "\t\t",gid
		print "determining which pairs or maps need to be compared"

	### get fits files from neighbors
	fits = [ (graceid, fit) for fit in fits ]
	for gid in graceids:
		if gid != graceid:
			fits += [ (gid, _) for _ in gracedb.files( gid ).json().keys() if _.endswith("fits") or _.endswith("fits.gz") ]

	fits.sort(key=lambda l: l[1])
	fits.sort(key=lambda l: l[0])

	### figure out which pairs have already happened
	compared = defaultdict( int )
	for gid in graceids:
		for log in gracedb.logs(gid).json()['log']:
			for ind, (gid1,fit1) in enumerate(fits):
				for (gid2,fit2) in fits[ind+1:]:
					compared[(gid1,fit1, gid2,fit2)] += "%s-%s,%s-%s : "%(gid1,fit1, gid2,fit2) in log['comment']

	### download files to local directory
	compare_cmds = []
	for ind,  (gid1, fit1) in enumerate(fits):
		for gid2, fit2 in fits[ind+1:]:
			if compared[(gid1,fit1, gid2,fit2)] == 0: ### no log messages found
				name1 = "%s-%s"%(gid1,fit1)
				if name1 not in names:
					if opts.verbose:
						print "\tdownloading %s from %s to %s"%(fit1, gid1, name1)
					file_obj = open(name1, "w")
					file_obj.write( gracedb.files( gid1, fit1).read() )
					file_obj.close()
					names.append( name1 )

				name2 = "%s-%s"%(gid2,fit2)
				if name2 not in names:
					if opts.verbose:
						print "\tdownloading %s from %s to %s"%(fit2, gid2, name2)
					file_obj = open(name2, "w")
					file_obj.write( gracedb.files( gid2, fit2).read() )
					file_obj.close()
					names.append( name2 )

				### add to compare command
				compare_cmds.append( compare_cmd + " %s,%s --graceid %s %s,%s --graceid %s "%(name1,name1, gid1, name2,name2, gid2) )

	### launch commands
	for compare_command in compare_cmds:
		if opts.verbose:
			print compare_command
		proc = subprocess.Popen(compare_command.split())
		proc.wait() # block!
		if opts.verbose:
			print "\tdone!"

#=================================================
### clean up files
if not opts.keep_maps:
	for name in names:
		if opts.verbose:
			print "removing", name
		os.remove( name )


