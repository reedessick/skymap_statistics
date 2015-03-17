usage = "grinch_executable.py [--options] graceid"
description = "manages the interface between GraceDB and the command line inputs for skymap quantification code"
author = "R. Essick (reed.essick@ligo.org)"

#=================================================

from ligo.gracedb.rest import GraceDB

import subprocess

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

gdburl = config.get('general', 'gdb_url')
tag_as_sky_loc = config.getbool('general', 'tag_as_sky_loc')

degrees = config.getbool('general', 'degrees')

credible_intervals = config.get('general', 'credible_interval').split()

#=================================================
### get neighbors from GraceDB
p_dt = config.getfloat('general', '+dt')
m_dt = config.getfloat('general', '-dt')

### command to get neighbors, then parse
neighbors = []

#=================================================
### need to query GraceDB to get fits files from events (including neighbors).

#=================================================
### need to define labels, graceid's, etc and pass them along to other executables

### only launch analyze_maps.py on the graceid specified as an argument

### launch compare_maps on all graceids, including neighbors

"""
don't upload information from within analyze and compare, but rather dump the statements into text files.
upload/tag those text files?

this will let us check whether we need to analyze particular files associated with a given event, or whether that has already been done?
"""
