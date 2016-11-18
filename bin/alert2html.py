#!/usr/bin/python
usage       = "alert2html.py [--options] config.ini"
description = "an lvalert_listen interface to generate and upload html summaries of FITS files"
author      = "reed.essick@ligo.org"

#-------------------------------------------------

import os
import sys

import json

import subprocess as sp

from ligo.gracedb.rest import GraceDb

from ConfigParser import SafeConfigParser

from optparse import OptionParser

#-------------------------------------------------

parser = OptionParser(usage=usage, description=description)

parser.add_option('-v', '--verbose', default=False, action='store_true')

parser.add_option('-g', '--graceid-and-filename', nargs=2, default=None, type='string',
    help='use for testing purposes (does not listen for alert thru stdin). eg: "G12345 bayestar.fits.gz' )

parser.add_option('', '--skip-gracedb-upload', default=False, action='store_true')

opts, args = parser.parse_args()

if len(args)!=1:
    raise ValueError('please supply exactly one input argument\n%s'%usage)
configname = args[0]

#-------------------------------------------------

### read in GraceID and Filename

if opts.graceid==None: ### listen for alert through stdin
    alert = sys.stdin.read()
    if opts.verbose:
        print( 'recieved alert : %s'%alert )
    alert = json.loads( alert )

    graceid = alert['uid'] ### pull out the graceid

    ###   must be an update alert                                    must be a FITS file
    if (alert['alert_type']=='update') and (alert['file'].endswith('.fits') or alert['file'].endswith('.fits.gz')):
        filename = alert['file'] ### pull out the filename

    else:
        print( 'ignoring...' )
        sys.exit(0)

else: ### graceid was supplied
    graceid, filename = opts.graceid_and_filename ### read in both graceid and filename      
        
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
json_nside = config.getint('general', 'json nside')

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
time_delay_alaph = config.get('time delay', 'alpha')

#------------------------

mollweide_levels     = config.get('mollweide', 'levels').split()
mollweide_alpha      = config.get('mollweide', 'alpha')
mollweide_linewidths = config.get('mollweide' 'linewidths')

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
marker_size      = conifg.get('marker', 'size')
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

### download FITS file
localname = os.path.join( outdir, filename )
if opts.verbose:
    print( 'downloading from %s : %s -> %s'%(graceDbURL, filename, localname) )

file_obj = open(localname, 'w')
file_obj.write( gracedb.files( graceid, filename ).read() )
file_obj.close()

#------------------------

### figure out which IFOs participated
ifos = gracedb.event( graceid ).json()['instruments'].split(',')

### format like I like them in this repo...
ifos = [ifo[0] for ifo in ifos] ### eg: H1 -> H

#-------------------------------------------------

### set up snglFITShtml command
snglFITScmd = [ 'snglFITShtml.py', localname,
                '--graceid',              graceid,
                '--graceDbURL',           graceDbURL,
                '--output-dir',           outdir,
                '--output-url',           outurl,
		'--json-nside',           json_nside,
                '--tag',                  tag,
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
    snglFITScmd += ['--conf', cnf]

### launch snglFITShtml
if opts.verbose:
    print( "    %s"%(" ".join(snglFITScmd)) )

proc = sp.Popen(snglFITScmd, stderr=sp.PIPE)
_, err = proc.communicate()
if proc.returncode:
    raise NotImplementedError('post a warning to GraceDb using error message from process? Or just print that here?')

#------------------------

### set up multFITShtml
multFITScmd = [ 'multFITShtml.py', 
                '--graceid',              graceid,
                '--graceDbURL',           graceDbURL,
                '--output-dir',           outdir,
                '--output-url',           outurl,
                '--tag',                  tag,
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

### search for other FITS files in outdir
multFITScmd += glob.glob( '%s/*.fits'%outdir ) + glob.glob( '%s/*.fits.gz'%outdir )

### launch multFITShtml
if opts.verbose:
    print( "    %s"%(" ".join(multFITScmd)) )

proc = sp.Popen(multFITScmd, stderr=sp.PIPE())
_, err = proc.communicate()
if proc.returncode:
    raise NotImplementedError('post a warning to GraceDb using error message from process? or just print that here?')
