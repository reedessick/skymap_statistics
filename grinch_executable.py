usage = "grinch_executable.py [--options] graceid"
description = "manages the interface between GraceDB and the command line inputs for skymap quantification code"
author = "R. Essick (reed.essick@ligo.org), A. Urban (aurban@uwm.edu)"

#=================================================

import json

import subprocess

from ligo.gracedb.rest import GraceDB

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

config = ConfigParser.SafeConfigParser()
config.read(opts.config)

### used to upload summary files and tag them to GraceDB
gdburl = config.get('general', 'gdb_url')
tag_as_sky_loc = config.getbool('general', 'tag_as_sky_loc')

### get universal options for all executables
universal_cmd = " ".join(" -c %.3f "%conf for conf in config.get('general', 'credible_interval') )
if config.getbool('general', 'degrees'):
	universal_cmd += " --degrees "

#=================================================
### get time for this GraceID
gevent = json.loads( gracedb.event(graceid).read() )

#=================================================
### get neighbors from GraceDB

### pull time windows from config file
p_dt = config.getfloat('general', '+dt')
m_dt = config.getfloat('general', '-dt')

### get time from GraceDB event
### query taken from A. Urban's raven/grace.py module
neighbors = [ json.load( _.read() ) for _ in list(gracedb.events( "%d..%d"%(t-m_dt, t+p_dt) )) ]

#=================================================
### need to query GraceDB to get fits files from events (including neighbors).

''' ### stolen from A. Urban's raven/grace.py module
def get_fits(gw_event):
    """ Downloads and unzips .fits file from gracedb into the 
        current working directory """
    os.system('gracedb download ' + gw_event.graceid + ' skymap.fits.gz')
'''

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

