#!/usr/bin/python
usage       = "alert2html.py [--options] config.ini"
description = "an lvalert_listen interface to generate and upload html summaries of FITS files"
author      = "reed.essick@ligo.org"

#-------------------------------------------------

import os
import sys
import glob

import json

import subprocess as sp

import fits2html

from ligo.gracedb.rest import GraceDb

from ConfigParser import SafeConfigParser

from optparse import OptionParser

#-------------------------------------------------

parser = OptionParser(usage=usage, description=description)

parser.add_option('-v', '--verbose', default=False, action='store_true')
parser.add_option('-V', '--Verbose', default=False, action='store_true')

parser.add_option('-g', '--graceid', default=None, type='string',
    help='download all fits files for this graceid and process them together. Forces multFITS to be tagged sky_loc \
(may produce redundant log messages if you are not careful)' )

parser.add_option('', '--graceid-and-filename', nargs=2, default=None, type='string',
    help='use for testing purposes (does not listen for alert thru stdin). eg: "G12345 bayestar.fits.gz". \
if --graceid is supplied, this will be ignored.' )

parser.add_option('', '--skip-gracedb-upload', default=False, action='store_true')

opts, args = parser.parse_args()

if len(args)!=1:
    raise ValueError('please supply exactly one input argument\n%s'%usage)
configname = args[0]

opts.verbose = opts.verbose or opts.Verbose

#-------------------------------------------------

### read in GraceID and Filename

if (opts.graceid==None) and (opts.graceid_and_filename==None): ### listen for alert through stdin

    alert = sys.stdin.read()
    if opts.Verbose:
        print( 'recieved alert : %s'%alert )
    alert = json.loads( alert )

    graceid = alert['uid'] ### pull out the graceid

    ###          must be an update alert                                    must be a FITS file
    if not ((alert['alert_type']=='update') and (alert['file'].endswith('.fits') or alert['file'].endswith('.fits.gz'))):
        if opts.Verbose:
            print( 'ignoring...' )
        sys.exit(0)

elif opts.graceid: ### graceid was supplied
    graceid = opts.graceid

else: ### --graceid was not supplied, but --graceid-and-filename was
    graceid = opts.graceid_and_filename[0] ### read in both graceid and filename      
        
#-------------------------------------------------

### read in Config file

if opts.verbose:
    print( 'reading config from : %s'%configname )
config = SafeConfigParser()
config.read(configname)

#------------------------

### extract some common parameters
graceDbURL = config.get('general', 'gracedb url')

### init connection to GraceDb
gracedb = GraceDb( graceDbURL )

#------------------------

### output placement and labeling
outdir = os.path.join( config.get('general', 'output dir'), graceid )
outurl = os.path.join( config.get('general', 'output url'), graceid )

tag = config.get('general', 'tag')
json_nside = config.get('general', 'json nside')

### make sure we have a working directory
if not os.path.exists(outdir):
    os.makedirs(outdir)

#------------------------

### output formatting
figtype = config.get('plotting', 'figtype')
dpi     = config.get('plotting', 'dpi')

transparent  = config.getboolean('plotting', 'transparent')
no_margticks = config.getboolean('plotting', 'no margticks')
color_map    = config.get('plotting', 'color map')

#------------------------

time_delay_color = config.get('time delay', 'color')
time_delay_alpha = config.get('time delay', 'alpha')

#------------------------

mollweide_levels     = config.get('mollweide', 'levels').split()
mollweide_alpha      = config.get('mollweide', 'alpha')
mollweide_linewidths = config.get('mollweide', 'linewidths')

#------------------------

continents_color = config.get('continents', 'color')
continents_alpha = config.get('continents', 'alpha')

#------------------------

line_of_sight_color = config.get('line of sight', 'color')

#------------------------

zenith_color = config.get('zenith', 'color')

#------------------------

marker           = config.get('marker', 'marker')
marker_color     = config.get('marker', 'color')
marker_alpha     = config.get('marker', 'alpha')
marker_size      = config.get('marker', 'size')
marker_edgewidth = config.get('marker', 'edgewidth')

#------------------------

dT_Nsamp   = config.get('dT marginals', 'Nsamp')
dT_nside   = config.get('dT marginals', 'nside')
dT_xlim_dB = config.get('dT marginals', 'xlim dB')

#------------------------

base  = config.get('stats', 'base')
confs = config.get('stats', 'conf').split()
areas = config.get('stats', 'area').split()

#-------------------------------------------------

### figure out which IFOs participated
ifos = gracedb.event( graceid ).json()['instruments'].split(',')

### format like I like them in this repo...
ifos = [ifo[0] for ifo in ifos] ### eg: H1 -> H

#-------------------------------------------------

### figure out which fits files we're analyzing
if opts.graceid: ### download all of them
    filenames = sorted( [ filename for filename in gracedb.files( graceid ).json().keys() if filename.endswith('.fits.gz') or filename.endswith('.fits') ] )

elif opts.graceid_and_filename: ### just this one
    filenames = [ opts.graceid_and_filename[1] ]

else: ### read from stdin
    filenames = [ alert['file'] ] ### pull out the filename

### search for other FITS files in outdir
localnames = sorted( glob.glob( '%s/*.fits'%outdir ) + glob.glob( '%s/*.fits.gz'%outdir ) )

if opts.Verbose:
    print "processing %d new FITS files"%len(filenames)
    for filename in filenames:
        print "    "+filename

    print "found %d existing FITS files"%len(localnames)
    for filename in localnames:
        print "    "+filename

#-------------------------------------------------
# BEGIN THE ANALYSIS
#-------------------------------------------------

### iterate over filenames, spawning a snglFITS job for each
for filename in filenames:

    if not opts.skip_gracedb_upload:
        gracedb.writeLog( graceid, message='started skymap summary for <a href="%s">%s</a>'%(os.path.join(graceDbURL, '..', 'events', graceid, 'files', filename), filename), tagname=fits2html.standard_tagname )

    #--------------------

    ### download FITS file
    localname = os.path.join( outdir, filename )
    if opts.verbose:
        print( 'downloading from %s : %s -> %s'%(graceDbURL, filename, localname) )

    file_obj = open(localname, 'w')
    file_obj.write( gracedb.files( graceid, filename ).read() )
    file_obj.close()

    #---------------------------------------------

    ### set up snglFITShtml command
    snglFITScmd = [ 'snglFITShtml.py', localname,
                '--graceDb-html-tagname', 'sky_loc', ### always tag this skyloc
                '--graceDbURL',           graceDbURL,
                '--output-dir',           outdir,
                '--output-url',           outurl,
		'--json-nside',           json_nside,
                '--figtype',              figtype,
                '--dpi',                  dpi,
                '--color-map',            color_map,
                '--time-delay-color',     time_delay_color,
                '--time-delay-alpha',     time_delay_alpha,
                '--mollweide-alpha',      mollweide_alpha,
                '--mollweide-linewidths', mollweide_linewidths,
                '--continents-color',     continents_color,
                '--continents-alpha',     continents_alpha,
                '--line-of-sight-color',  line_of_sight_color,
                '--zenith-color',         zenith_color,
                '--marker',               marker,
                '--marker-color',         marker_color,
                '--marker-alpha',         marker_alpha,
                '--marker-size',          marker_size,
                '--marker-edgewidth',     marker_edgewidth,
                '--dT-Nsamp',             dT_Nsamp,
                '--dT-nside',             dT_nside,
                '--dT-xlim-dB',           dT_xlim_dB,
                '--base',                 base,
              ]

    if tag:
        snglFITScmd += ['--tag', tag]

    if opts.skip_gracedb_upload:
        snglFITScmd.append( '--skip-gracedb-upload' )
    else:
        snglFITScmd += ['--graceid', graceid ]

    if opts.verbose:
        snglFITScmd.append( '--verbose' )

    if transparent:
        snglFITScmd.append( '--transparent' )

    if no_margticks:
        snglFITScmd.append( '--no-margticks' )

    for ifo in ifos:
        snglFITScmd += ['-i', ifo]

    for level in mollweide_levels:
        snglFITScmd += ['--mollweide-levels', level]

    for conf in confs:
        snglFITScmd += ['--conf', conf]

    ### launch snglFITShtml
    if opts.verbose:
        print( "    %s"%(" ".join(snglFITScmd)) )

    proc = sp.Popen(snglFITScmd, stderr=sp.PIPE)
    _, err = proc.communicate()
    if proc.returncode:
        print err
        if opts.skip_gracedb_upload:
            raise RuntimeError('snglFITS returncode=%d\n%s'%(proc.returncode, " ".join(snglFITScmd)))

        else:
            print( 'WARNING: snglFITS returncode=%d'%proc.returncode )
            gracedb.writeLog( graceid, 'WARNING: snglFITS failed for <a href="%s">%s</a>'%(os.path.join(graceDbURL, '..', 'events', graceid, 'files', filename), filename), tagname=fits2html.standard_tagname )

    elif not opts.skip_gracedb_upload:
        gracedb.writeLog( graceid, message='finished skymap summary for <a href="%s">%s</a>'%(os.path.join(graceDbURL, '..', 'events', graceid, 'files', filename), filename), tagname=fits2html.standard_tagname )


#------------------------

### set up multFITShtml

all_localnames = sorted(set(localnames + [os.path.join( outdir, filename ) for filename in filenames]))

if not opts.skip_gracedb_upload:
    filename = all_localnames[0]
    localname = os.path.basename( filename )
    hrefs = '<a href="%s">%s</a>'%(os.path.join(graceDbURL, '..', 'events', graceid, 'files', localname), localname)

    for filename in all_localnames[1:]:
        localname = os.path.basename( filename )
        hrefs += ', <a href="%s">%s</a>'%(os.path.join(graceDbURL, '..', 'events', graceid, 'files', localname), localname)

    message = 'started skymap comparison for %s'%hrefs
    gracedb.writeLog( graceid, message=message, tagname=fits2html.standard_tagname )

multFITScmd = [ 'multFITShtml.py', 
                '--graceDbURL',           graceDbURL,
                '--output-dir',           outdir,
                '--output-url',           outurl,
                '--figtype',              figtype,
                '--dpi',                  dpi,
                '--color-map',            color_map,
                '--time-delay-color',     time_delay_color,
                '--time-delay-alpha',     time_delay_alpha,
                '--mollweide-alpha',      mollweide_alpha,
                '--mollweide-linewidths', mollweide_linewidths,
                '--continents-color',     continents_color,
                '--continents-alpha',     continents_alpha,
                '--line-of-sight-color',  line_of_sight_color,
                '--zenith-color',         zenith_color,
                '--marker',               marker,
                '--marker-color',         marker_color,
                '--marker-alpha',         marker_alpha,
                '--marker-size',          marker_size,
                '--marker-edgewidth',     marker_edgewidth,
                '--dT-Nsamp',             dT_Nsamp,
                '--dT-nside',             dT_nside,
                '--dT-xlim-dB',           dT_xlim_dB,
                '--base',                 base,
              ]

multFITScmd += all_localnames

if tag:
    multFITScmd += ['--tag', tag]

if opts.skip_gracedb_upload:
    multFITScmd.append( '--skip-gracedb-upload' )
else:
    multFITScmd += ['--graceid', graceid]

if opts.verbose:
    multFITScmd.append( '--verbose' )

if transparent:
    multFITScmd.append( '--transparent' )

if no_margticks:
    multFITScmd.append( '--no-margticks' )

for ifo in ifos:
    multFITScmd += ['-i', ifo]

for level in mollweide_levels:
    multFITScmd += ['--mollweide-levels', level]

for conf in confs:
    multFITScmd += ['--conf', conf]

for area in areas:
    multFITScmd += ['--area', area]

if len(localnames)==0 or opts.graceid: ### no previous FITS files, so this is the first skymap (FIXME: subject to a race condition, but should be pretty unlikely)
                                       ### or, we supplied --graceid so we want to force this action
    multFITScmd += ['--graceDb-html-tagname', 'sky_loc']

### launch multFITShtml
if opts.verbose:
    print( "    %s"%(" ".join(multFITScmd)) )

proc = sp.Popen(multFITScmd, stderr=sp.PIPE)
_, err = proc.communicate()
if proc.returncode:
    print err
    if opts.skip_gracedb_upload:
        raise NotImplementedError('WARNING: multFITS returncode=%d\n%s'%(proc.returncode, " ".join(multFITScmd)))
    else:
        print( 'WARNING: multFITS returncode=%d'%proc.returncode )
        gracedb.writeLog( graceid, 'WARNING: multFITS failed for <a href="%s">%s</a>'%(os.path.join(graceDbURL, '..', 'events', graceid, 'files', filename), filename), tagname=fits2html.standard_tagname )

#------------------------

if not opts.skip_gracedb_upload:
    message = 'finished skymap comparison for %s'%hrefs
    gracedb.writeLog( graceid, message=message, tagname=fits2html.standard_tagname )
