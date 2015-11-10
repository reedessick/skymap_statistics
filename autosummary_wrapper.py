#!/usr/bin/python
usage = "launch_autosummary.py [--options]"
description = "a wrapper around autosummary.py to set it up and run under lvalert_listen"
auther = "Reed Essick (reed.essick@ligo.org)"

import sys
import json

import subprocess as sp

from ConfigParser import SafeConfigParser
from optparse import OptionParser

#=================================================

parser = OptionParser(usage=usage, description=description)

parser.add_option("-v", "--verbose", default=False, action="store_true")

parser.add_option("-c", "--config", default="./config.ini", type="string")

opts, args = parser.parse_args()

#=================================================

### read lvalert 
alert = sys.stdin.read()

if opts.verbose:
    print "alert received:\n%s"%alert

alert = json.loads(alert)

### determine if we need to react (only when there is a new FITS file)
if (alert['alert_type'] == 'update') and alert['filename'].strip(".gz").endswidth(".fits"):  ### check for new FITS file 

    ### configure command
    if opts.verbose:
        print "new FITS file : %s -> %s"%(alert['uid'], alert['filename'])
        print "reading config : %s"%(opts.config)
    config = SafeConfigParser()
    config.read( opts.config )

    options = dict( config.items( "general" ) )
  
    event_type = None
    if config.has_section( event_type ):
        if opts.verbose:
            print "\tloading extra instructions from section : %s"%event_type
        options.update( dict( config.items( event_type ) ) )
    elif verbose:
        print "\tno section found for event_type : %s"%event_type



    cmd = "autosummary.py %s
    ### launch

