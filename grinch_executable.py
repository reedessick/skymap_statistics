usage = "grinch_executable.py [--options] graceid"
description = "manages the interface between GraceDB and the command line inputs for skymap quantification code"
author = "R. Essick (reed.essick@ligo.org), A. Urban (aurban@uwm.edu)"

#=================================================

from math import ceil, floor
import json

import subprocess

from ligo.gracedb.rest import GraceDb

import ConfigParser
from optparse import OptionParser

#=================================================

parser = OptionParser(usage=usage, description=description)

parser.add_option("-v", "--verbose", default=False, action="store_true")

parser.add_option("-c", "--config", default="grinch_config.ini", type="string")

opts, args = parser.parse_args()

if len(args)!= 1:
	raise ValueError("please supply exactly one GraceID as an argument!")
graceid = args[0]

#=================================================

if opts.verbose: 
	print "reading config from", opts.config

config = ConfigParser.SafeConfigParser()
config.read(opts.config)

### used to upload summary files and tag them to GraceDB
gdburl = config.get('general', 'gdb_url')
tag_as_sky_loc = config.getboolean('general', 'tag_as_sky_loc')

### instantiate gracedb interface
if opts.verbose: 
	print "instantiating connection to GraceDb :", gdburl
gracedb = gracedb = GraceDb( gdburl )

#=================================================

### get universal options for all executables
universal_cmd = " ".join(" -c %s "%(conf) for conf in config.get('general', 'credible_interval').split() )
if config.getboolean('general', 'degrees'):
	universal_cmd += " --degrees "

#=================================================
### grab gracedb events and their files!
gdb_entries = {}

#========================
### get time for this GraceID
if opts.verbose: 
	print "reading data for :", graceid

gdb_entries[graceid] = gracedb.event(graceid).json()

#========================
### get neighbors from GraceDB

t = float(gdb_entries[graceid]['gpstime'])

### pull time windows from config file
p_dt = config.getfloat('general', '+dt')
m_dt = config.getfloat('general', '-dt')

### get time from GraceDB event
if opts.verbose: 
	print "searching for neighbors within [%.5f, %.5f]"%(t+m_dt, t+p_dt)
gdb_entries.update( dict( (gdb_entry['graceid'], gdb_entry) for gdb_entry in gracedb.events( "%d..%d"%(floor(t+m_dt), ceil(t+p_dt)) ) if gdb_entry['graceid']!=graceid ) )

### downselect this based on event type?

if opts.verbose:
	print "\tfound %d neighbors"%(len(gdb_entries)-1)
	for gid in gdb_entries.keys():
		if gid != graceid:
			print "\t\t", gid 

#=================================================
### need to query GraceDB to get fits files from events (including neighbors).

if opts.verbose:
	print "querying for files"

todo = {}
for gid in gdb_entries.keys():
	files = gracedb.files(gid).json()
	fits = [ _ for _ in files if _.endswith("fits") or _.endswith("fits.gz") ]
	todo.update( dict( (gid, fit) for fit in fits ) )


### need to define a standard filename to determine whether we have analyzed a given fits file!

#=================================================
### need to define labels, graceid's, etc and pass them along to other executables
### with current work flow, how do we prevent multiple entries/redundant processing?
### 	how do we control that when we don't know what each executable does when launching (desireable for code simplicity?)

### only launch analyze_maps.py on the graceid specified as an argument
### launch compare_maps on all graceids, including neighbors
### launch these in parallel?

# iterate over executables
for section, cmd in config.items('executables'):
	if not config.has_section(section): ### check that section exists!
		raise ValueError("no section \"%s\" in %s"%(section, opts.config) )

	### extract options and build cmd
	cmd = "%s %s %s"%(cmd, universal_cmd, " ".join(" --%s "%option for option, boolean in config.items(section) if bool(boolean)))

	if opts.verbose:
		print cmd

	### add "dont_wait" functionality?
#	proc = subprocess.Popen(cmd.split())
#	proc.wait()

"""
don't upload information from within analyze and compare, but rather dump the statements into text files.
upload/tag those text files?

this will let us check whether we need to analyze particular files associated with a given event, or whether that has already been done?
"""

