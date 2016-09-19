description = "a module that houses classes that write html structures for fits2html.py and friends"
author      = "reed.essick@ligo.org"

#-------------------------------------------------

import os
import json

import stats
import triangulate

from plotting import mollweide as mw
from plotting import cartesian as ct
from plotting import colors

plt = mw.plt

import detector_cache
import antenna

import numpy as np
import healpy as hp

from lal.gpstime import tconvert

from ligo.gracedb.rest import GraceDb

#-------------------------------------------------

class Figure(object):
    '''
    a thin wrapper around figure objects that knows how to upload and save them
    '''

    def __init__(self, fig, output_dir, output_url, graceid=None, graceDbURL='https://gracedb.ligo.org/api/'):
        self.fig = fig

        self.output_dir = output_dir
        self.output_url = output_url

        self.graceid    = graceid
        self.graceDbURL = graceDbURL

    def saveAndUpload(self, figname, message='', tagname=['skymapAutosummary']):
        figname = os.path.join(output_dir, figname)
        self.fig.savefig( figname )
        if self.graceid!=None:
            gdb = GraceDb(self.graceDbURL)
            httpResponse = gdb.writeLog( self.graceid, message=message, filename=figname, tagname=tagname )
            ### may want to check httpResponse for errors...

        return "%s/%s"%(self.output_url, figname) 

class Json(object):
    '''
    a thin wrapper around a json object that knows how to upload and save them
    '''

    def __init__(self, obj, output_dir, output_url, graceid=None, graceDbURL='https://gracedb.ligo.org/api/'):
        self.obj = obj

        self.output_dir = output_dir
        self.output_url = output_url

        self.graceid    = graceid
        self.graceDbURL = graceDbURL

    def saveAndUpload(self, filename, message='', tagname=['skymapAutosummary']):
        filename = os.path.join(output_dir, filename)
        file_obj = open( filename, "w" )
        json.dump( self.obj, file_obj )
        file_obj.close()
        if self.graceid!=None:
            gdb = GraceDb(self.graceDbURL)
            httpResponse = gdb.writeLog( self.graceid, message=message, filename=filename, tagname=tagname )
            ### may want to check httpResponse for errors...

        return "%s/%s"%(self.output_url, figname)

#-------------------------------------------------

class snglFITS(object):
    '''
    a class that houses data and renders html, json for info about a single FITS file (ie: no comparisons)
    '''

    def __init__( self, 
                  fitsname, 
                  ### general options about output
                  output_dir = '.',
                  output_url = './', ### ignored if graceid is supplied, otherwise used to build URLs in html document
                  tag        = '',
                  figtype    = "png",
                  dpi        = 500,
                  graceid    = None, ### if supplied, upload files and reference them in the html document
                  graceDbURL = 'https://gracedb.ligo.org/api/',
                  ### general options about annotation and which plots to build
                  ifos = [],
                  ### general options about colors, shading, and labeling
                  color_map   = "OrRd",
                  transparent = False,
                  no_yticks   = False,
                  ### options about mollweide projections
                  mollweide_levels    = [0.1, 0.5, 0.9],
                  mollweide_alpha     = 1.0, 
                  mollweide_linewidth = 1.0,
                  time_delay_color    = 'k',
                  time_delay_alpha    = 1.0, 
                  line_of_sight_color = 'k',
                  zenith_color        = 'k',
                  marker              = 'o', 
                  marker_color        = 'k', 
                  marker_alpha        = 1.0, 
                  marker_size         = 4, 
                  marker_edgewidth    = 1, 
                  continents          = True,
                  continents_color    = 'k',
                  continents_alpha    = 0.5,
                  ### plotting options for dT marginals
                  dT_Nsamp     = 1001,
                  dT_nside     = None,
                  dT_xlim_dB   = -20,
                  ### options for computing statistics
                  base = 2.0,
                  conf = np.linspace(0,1,51), 
                ):

        ### general things about FITS file
        self.fitsname = fitsname
        self.gps      = gps

        self.readFITS() ### read in FITS and store local copies

        ### which IFOs are important
        for ifo in ifos:
            assert detector_cache.detectors.has_key(ifo), "ifo=%s is not understood!"%ifo
        self.ifos = sorted(ifos)
        self.ifo_pairs = []
        for ind, ifo1 in enumerate(self.ifos):
            for ifo2 in self.ifos[ind+1:]:
                self.ifo_pairs.append( (ifo1, ifo2) )
        
        ### output formatting
        self.output_dir = output_dir
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        if graceid!=None:
            output_url = graceDbURL+'../events/%s/files/'%(graceid)
        else:
            output_url = output_url

        self.tag = tag

        self.figtype = figtype
        self.dpi     = dpi

        self.graceid    = graceid
        self.gracedbURL = graceDbURL

        ### general color schemes
        self.color_map   = color_map
        self.transparent = transparent
        self.no_yticks   = no_yticks

        ### options for mollweide projections
        self.mollweide_levels    = mollweide_levels
        self.mollweide_alpha     = mollweide_alpha 
        self.mollweide_linewidth = mollweide_linewidth

        self.time_delay_color = time_delay_color
        self.time_delay_alpha = time_delay_alpha

        self.line_of_sight_color = line_of_sight_color
        self.zenith_color        = zenith_color

        self.marker           = marker
        self.marker_color     = marker_color
        self.marker_alpha     = marker_alpha
        self.marker_size      = marker_size
        self.marker_edgewidth = marker_edgewidth
        
        self.continents       = continents
        self.continents_color = continents_color
        self.continents_alpha = continents_alpha

        ### options for time-delay marginals
        self.dT_Nsamp     = dT_Nsamp
        self.dT_nside     = dT_nside
        self.dT_xlim_dB   = dT_xlim_dB
        
        ### options for statistics
        self.base = base 
        self.conf = conf

        ### local references for plotting
        self.figind = 0

    def readFITS(self, verbose=False):
        '''
        reads in the FITS file and sets up local copies
        '''
        ### load in map
        post, header = hp.read_map( self.fitsname, h=True, verbose=verbose )
        header = dict(header)

        ### ensure we are in RING ordering
        if h['ORDERING']=='NEST':
            post = hp.reorder( post, h['ORDERING'], r2n=1 )

        ### extract gps time
        self.gps = tconvert(h['DATE-OBS'])

        ### set up references to maps in C and E coordinates
        if verbose:
            print "  setting up local copies in C and E coordinates"
        coord = h['COORDSYS']
        if coord == 'C':
            self.postC = post[:]
            self.postE = triangulate.rotateMapC2E( postC, self.gps )
        elif coord == 'E':
            self.postE = post[:]
            self.postC = triangulate.rotateMapE2C( postE, self.gps )
        else:
            raise ValueError('COORDSYS=%s not understood!'%coord)

        ### set up meta-data about the map
        if verbose:
            print "  setting up local references to angles"
        npix = len(post)
        self.nside = hp.npix2nside( npix )
        self.theta, self.phi = hp.pix2ang( self.nside, np.arange(npix) )

        ### compute basic information statistics about map
        self.entropy     = stats.entropy( post, base=self.base )
        self.information = stats.information( post, base=self.base )

    def make_mollweide(self, verbose=False):
        """
        make mollweide projections
        """
        if verbose:
            print "building mollweide projections"
        self.mollweide = dict()

        for projection, coord, post in [('astro hour mollweide', 'C', self.postC), ('mollweide', 'E', self.postE)]:

            ### generate figure
            fig, ax = mw.gen_fig_ax( self.figind, projection=projection )
            fig = Figure( fig, self.output_dir, self.output_url, graceid=self.graceid, graceDbURL=self.graceDbURL )
            self.figind += 1

            mw.heatmap( post, ax, color_map=opts.color_map )
            mw.annotate( ax,
                         continents       = self.continents and (coord=='E'),
                         continents_color = self.continents_color,
                         continents_alpha = self.continents_alpha,
                       )

            if self.transparent:
                fig.fig.patch.set_alpha(0.)
                ax.patch.set_alpha(0.)
                ax.set_alpha(0.)

            ### save just the heatmap
            figname = "heatmap%s%s.%s"%(coord, self.tag, self.figtype)
            if verbose:
                print "  "+figname
            self.mollweide[coord] = fig.saveAndUpload( figname )

            ### annotate with fancy crap
            mapIND = post.argmax()
            mapY = self.theta[mapIND]
            mapX = self.phi[mapIND]
            if coord == 'C':
                mapY = 0.5*np.pi - mapY ### convert from Theta->Dec

            mw.annoate( ax,
                        projection          = projection,
                        line_of_sight       = mw.gen_line_of_sight( self.ifo_pairs, coord=coord, gps=self.gps ),
                        line_of_sight_color = self.line_of_sight_color,
                        zenith              = mw.gen_zenith( self.ifo, coord=coord, gps=self.gps ),
                        zenith_color        = self.zenith_color,
                        time_delay          = mw.gen_time_delay( [(mapY, mapX)], self.ifo_pairs, coord=coord, gps=self.gps, degrees=False ),
                        time_delay_color    = self.time_delay_color,
                        time_delay_alpha    = self.time_delay_alpha,
                        marker_Dec_RA       = mw.gen_marker_Dec_RA( [(mapY, mapX)], coord=coord, gps=self.gps, degrees=False ),
                        marker              = self.marker,
                        marker_color        = self.marker_color,
                        marker_size         = self.marker_size,
                        marker_edgewidth    = self.marker_edgewidth,
                        marker_alpha        = self.marker_alpha,
                      )

            ### save heatmap + fancy crap
            figname = "heatmap%s-annotated%s.%s"%(coord, self.tag, self.figtype)
            if verbose:
                print "  "+figname
            self.mollweide[coord+" ann"] = fig.saveAndUpload( figname )

            ### add antenna patterns as contours
            genColor = colors.getColors()
            for ind, ifo in enumerate(self.ifo):

                Fp, Fx = detector_cache.detectors[ifo].antenna_patterns( self.theta, self.phi, 0.0 )
                ant = Fp**2 + Fx**2
                ant /= np.sum(ant)
                if coord == 'C':
                    ant = triangulate.rotateMapE2C( ant, self.gps )

                color = genColor.next()
                mw.contour( ant,
                            ax,
                            colors     = color,
                            levels     = self.mollweide_levels,
                            alpha      = self.mollweide_alpha,
                            linewidths = self.mollweide_linewidths 
                          )
                fig.text(0.01, 0.99-0.05*ind, ifo, color=color, ha='left', va='top')

            ### save heatmap + fancy crap + antenna pattern contours
            figname = "heatmap%s-antennas%s.%s"%(coord, self.tag, self.figtype)
            if verbose:
                print "  "+figname
            self.mollweide[coord+" ant"] = fig.saveAndUpload( figname )

            ### done with this figure
            plt.close( fig.fig )
            del fig

    def make_dT(self, verbose=False):
        '''
        make time-delay marginal plots and statistics
        '''
        if verbose:
            print "building time-delay marginals"
        self.dT = dict()

        for ifo1, ifo2 in self.ifo_pairs:
            ifos = "".join([ifo1, ifo2])
            if verbose:
                print "  %s - %s"%(ifo1, ifo2)

            sampDt = ct.get_sampDt( ifos, Nsamp=opts.dT_Nsamp )
            maxDt = sampDt[-1]

            fig, ax = ct.genDT_fig_ax( self.figind )
            fig = Figure( fig, self.output_dir, self.output_url, graceid=self.graceid, graceDbURL=self.graceDbURL )
            self.figind += 1

            ax.set_xlim(xmin=maxDt*1e3, xmax=-maxDt*1e3) ### we work in ms here...

            if self.dT_nside:
                kde = ct.post2marg( stats.resample(postE, self.dT_nside), ifos, sampDt, coord='E' )
            else:
                kde = ct.post2marg( postE, ifos, sampDt, coord='E' )

            ### compute statistics of the marginal
            raise NotImplementedError('compute statistics for dT marginals: entropy, information (in base)')

            ### plot
            ct.plot( ax, sampDt, kde, xlim_dB=self.dT_xlim_dB )

            ### decorate
            ax.set_xlabel(r'$\Delta t_{%s}\ [\mathrm{ms}]$'%(ifos))
            ax.set_ylabel(r'$p(\Delta t_{%s}|\mathrm{data})$'%(ifos))

            if opts.no_yticks:
                ax.set_yticklabels([])

            if self.transparent:
                fig.fig.patch.set_alpha(0.)
                ax.patch.set_alpha(0.)
                ax.set_alpha(0.)

            ### save just dT marginals
            figname = "dT_%s%s.%s"%(ifos, opts.tag, opts.figtype)
            if verbose:
                print "  "+figname
            self.dT['dT '+ifos] = fig.saveAndUpload( figname )

            ### annotate the plot
            ct.annotate( ax,
                         [ hp.pix2ang( nside, np.argmax(postE) ) ],
                         ifos,
                         maxDt,
                         coord   = 'E',
                         gps     = opts.gps,
                         color   = opts.time_delay_color,
                         alpha   = opts.time_delay_alpha,
                         degrees = opts.False,
                       )

            ### save annotated dT marginals
            figname = "dT_%s-annotated%s.%s"%(ifos, opts.tag, opts.figtype)
            if opts.verbose:
                print figname
            self.dT['dT ann '+ifos] = fig.saveAndUpload( figname )

            plt.close(fig.fig)
            del fig

    def make_los(self, verbose=False):
        '''
        make line-of-sight cartesian projections and statistics
        '''
        if verbose:
            print "building line-of-sight cartesian projections"
        self.los = dict() 

        for ifo1, ifo2 in self.ifo_pairs:
            if verbose:
                print "  %s - %s"%(ifo1, ifo2)

            t, p = triangulate.line_of_sight( ifo1, ifo2, coord='E' )

            fig, ax, rproj, tproj = ct.genHist_fig_ax( self.figind, figwidth=opts.figwidth, figheight=opts.figheight )
            fig = Figure( fig, self.output_dir, self.output_url, graceid=self.graceid, graceDbURL=self.graceDbURL )
            self.figind += 1

            ### rotate
            rtheta, rphi = triangulate.rotate2pole( self.theta, self.phi, t, p )

            Nbins = max(100, int(npix**0.5/5))

            ### compute mutual info
            raise NotImplementedError('compute mutual information distance')
#            mi, entj = triangulate.compute_mi( rtheta, rphi, Nbins, weights=self.postE )

            ### plot
            ct.histogram2d( rtheta, 
                            rphi, 
                            ax, 
                            rproj, 
                            tproj, 
                            Nbins   = Nbins, 
                            weights = self.postE, 
                            contour = self.contour, 
                            color   = colors.genColor().next(), ### always get the first color!
                            cmap    = self.color_map 
                          )

            ### save
            figname = "los-%s-%s%s.%s"%(ifo1, ifo2, self.tag, self.figtype)
            if verbose:
                print "  "+figname
            self.los['%s%s'%(ifo1,ifo2)] = fig.saveAndUpload( figname )

            plt.close( fig.fig )
            del fig

    def make_json(self, verbose=False):
        '''
        write map in C coords to json file
        '''
        raise NotImplementedError

    def make_cumulative_json(self, verbose=False):
        '''
        write cumulative map in C coords to json file
        '''
        raise NotImplementedError

    def make_confidence_regions(self, verbose=False):
        '''
        compute confidence regions
        '''
        raise NotImplementedError

    def make_antenna_patterns(self, verbose=False):
        '''
        compute antenna pattern statistics
        '''
        raise NotImplementedError


        '''
if opts.verbose:
    print "computing statistics"

data['nside'] = nside

### information measures
if opts.verbose:
    print "  information measures"
data['H'] = stats.entropy( postC, base=opts.base )
data['I'] = stats.information( postC, base=opts.base )

### size of confidence regions
if opts.verbose:
    print "  confidence regions"
pixarea = hp.nside2pixarea( nside )
cr = stats.credible_region(postC, opts.conf)
cr_size = [len(_)*pixarea for _ in cr]
data['CR size'] = zip( opts.conf, cr_size )

fig, ax = ct.genCR_fig_ax( figind, figwidth=opts.figwidth, figheight=opts.figheight )
figind += 1

ct.plot( ax, opts.conf, cr_size )

ax.set_xlabel('confidence')
ax.set_ylabel('confidence region [deg$^2$]')

figname = "%s/crSize%s.%s"%(outdir, opts.tag, opts.figtype)
if opts.verbose:
    print "    %s"%(figname)
fig.savefig( figname )
data['CR fig'] = figname

### antenna patterns
if opts.verbose:
    print "  antenna patterns"

mapind = np.argmax(postE)
for ifo in opts.ifo:
    if opts.verbose:
        print "    %s"%ifo

    Fp, Fx = detector_cache.detectors[ifo].antenna_patterns( theta, phi, 0.0 ) ### in Geographic coordinates

    ### MAP values
    mapind = np.argmax(postE)
    data['%s ant map'%ifo] = Fp[mapind]**2 + Fx[mapind]**2

    ### averaged over posterior
    data['%s ant ave'%ifo] = np.sum( postE*(Fp**2 + Fx**2) )

#-------------------------------------------------

### if we have posterior_samples.dat, generate postviz (only relevant for LIB and LALInf)

if opts.posterior_samples:
    raise NotImplementedError('parser posterior_samples and generate postviz page')
    ### use opts.posterior_samples_var to determine headers
    ### use opts.posterior_samples_decrement to get newNSIDE = nside/opts.posterior_samples_decrement, and use newNSIDE when (greedy) binning samples
'''





























    def toSTR(self):
        htmlSTR = """<head>
</head>
<body>
<p>%s<br>nside=%d<br>H=%.3f deg2<br>I=%.3f deg2</p>
"""%(fitsname, self['nside'], self['H'], self['I'])

        ### mollweide projections
        htmlSTR += "\n<hr>"
        for coord in "C E".split():
             htmlSTR += """
<img src=\"%s\"><img src=\"%s\"><img src=\"%s\"><br>"""%(self['mw %s'%coord][0], self['mw %s ann'%coord][0], self['mw %s ant'%coord][0])

        ### zenith frames
        htmlSTR += "\n<hr>"
        for ifo in self.ifos:
            htmlSTR += """
<img src=\"%s\">mi=%.3f, Hj=%.3f<br>"""%(self['zen %s'%ifo][0], self['zen %s mi'%ifo], self['zen %s Hj'])

        ### line-of-sight frames
        htmlSTR += "\n<hr>"
        for ind, ifo1 in enumerate(self.ifos):
            for ifo2 in self.ifos[ind+1:]:
                htmlSTR += """
<img src=\"%s\">mi=%.3f, Hj=%.3f<br>"""%(self['los %s%s'%(ifo1, ifo2)][0], self['los %s%s mi'%(ifo1,ifo2)], self['los %s%s Hj'%(ifo1,ifo2)])

        ### confidence regions
        htmlSTR += "\n<hr>"
        htmlSTR += "\n<img src=\"%s\">"%(self['CR fig'][0])
        for conf, size in self['CR size']:
            htmlSTR += "\n<p>conf=%.3f  size=%.3f deg2</p>"%(conf, size)

        ### antenna patterns
        htmlSTR += "\n<hr>"
        for ifo in self.ifos:
            htmlSTR += "<p>%s Fp^2+Fx^2 at MAP=%.3f, average=%.3f</p>"%(ifo, self['%s ant map'%ifo], self['%s ant ave'%ifo])

        htmlSTR += "\n</body>"

        return htmlSTR

    def write(self, filename):
        file_obj = open(filename, "w")
        file_obj.write( self.toSTR() )
        file_obj.close()

#---------------------------------------------------------------------------------------------------

class multFITS(object):
    '''
    a class that houses data and renders html, json for comparison info about multiple FITS files
    '''

    def __init__(self, *args):
        self.fitsnames = []
        super(multFITS, self).__init__(*args)

    def write(self, filename):
        raise NotImplementedError('should write an html document into filename=%s'%filename)
