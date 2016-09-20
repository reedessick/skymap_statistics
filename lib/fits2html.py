description = "a module that houses classes that write html structures for fits2html.py and friends"
author      = "reed.essick@ligo.org"

#-------------------------------------------------

from html import HTML
import getpass

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
        filename = os.path.join(self.output_dir, figname)
        self.fig.savefig( filename )
        if self.graceid!=None:
            gdb = GraceDb(self.graceDbURL)
            httpResponse = gdb.writeLog( self.graceid, message=message, filename=filename, tagname=tagname )
            ### may want to check httpResponse for errors...

        return os.path.join(self.output_url, figname) 

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
        fileName = os.path.join(self.output_dir, filename)
        file_obj = open( fileName, "w" )
        json.dump( self.obj, file_obj )
        file_obj.close()
        if self.graceid!=None:
            gdb = GraceDb(self.graceDbURL)
            httpResponse = gdb.writeLog( self.graceid, message=message, filename=fileName, tagname=tagname )
            ### may want to check httpResponse for errors...

        return os.path.join(self.output_url, filename)

#-------------------------------------------------

class snglFITS(object):
    '''
    a class that houses data and renders html, json for info about a single FITS file (ie: no comparisons)

    this class also encapsulates how the html pages are structured. This should be mirrored in multFITS as much as possible.
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
                  ### options for json reference files
                  json_nside = 128,
                  ### general options about annotation and which plots to build
                  ifos = [],
                  ### general options about colors, shading, and labeling
                  color_map   = "OrRd",
                  transparent = False,
                  no_yticks   = False,
                  ### options about mollweide projections
                  mollweide_levels     = [0.1, 0.5, 0.9],
                  mollweide_alpha      = 1.0, 
                  mollweide_linewidths = 1.0,
                  time_delay_color     = 'k',
                  time_delay_alpha     = 1.0, 
                  line_of_sight_color  = 'k',
                  zenith_color         = 'k',
                  marker               = 'o', 
                  marker_color         = 'k', 
                  marker_alpha         = 1.0, 
                  marker_size          = 4, 
                  marker_edgewidth     = 1, 
                  continents_color     = 'k',
                  continents_alpha     = 0.5,
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
        self.label    = os.path.basename(fitsname).split('.')[0]

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
            self.output_url = graceDbURL+'../events/%s/files/'%(graceid)
        else:
            self.output_url = output_url

        self.tag = tag

        self.figtype = figtype
        self.dpi     = dpi

        self.graceid    = graceid
        self.graceDbURL = graceDbURL

        ### general color schemes
        self.color_map   = color_map
        self.transparent = transparent
        self.no_yticks   = no_yticks

        ### options for mollweide projections
        self.mollweide_levels     = mollweide_levels
        self.mollweide_alpha      = mollweide_alpha 
        self.mollweide_linewidths = mollweide_linewidths

        self.time_delay_color = time_delay_color
        self.time_delay_alpha = time_delay_alpha

        self.line_of_sight_color = line_of_sight_color
        self.zenith_color        = zenith_color

        self.marker           = marker
        self.marker_color     = marker_color
        self.marker_alpha     = marker_alpha
        self.marker_size      = marker_size
        self.marker_edgewidth = marker_edgewidth
        
        self.continents_color = continents_color
        self.continents_alpha = continents_alpha

        ### options for time-delay marginals
        self.dT_Nsamp     = dT_Nsamp
        self.dT_nside     = dT_nside
        self.dT_xlim_dB   = dT_xlim_dB
        
        ### options for statistics
        self.base = base 
        self.conf = conf

        ### options for json reference files
        self.json_nside = json_nside

        ### local references for plotting
        self.figind = 0

    def readFITS(self, verbose=False):
        '''
        reads in the FITS file and sets up local copies
        '''
        if verbose:
            print "reading : "+self.fitsname
        ### load in map
        post, header = hp.read_map( self.fitsname, h=True, verbose=verbose )
        header = dict(header)

        ### ensure we are in RING ordering
        if header['ORDERING']=='NEST':
            post = hp.reorder( post, 'NEST', r2n=1 )

        ### extract gps time
        self.gps = tconvert(header['DATE-OBS'])

        ### set up references to maps in C and E coordinates
        if verbose:
            print "  setting up local copies in C and E coordinates"
        coord = header['COORDSYS']
        if coord == 'C':
            self.postC = post[:]
            self.postE = triangulate.rotateMapC2E( self.postC, self.gps )
        elif coord == 'E':
            self.postE = post[:]
            self.postC = triangulate.rotateMapE2C( self.postE, self.gps )
        else:
            raise ValueError('COORDSYS=%s not understood!'%coord)

        ### set up meta-data about the map
        if verbose:
            print "  setting up local references to angles"
        self.npix = len(post)
        self.nside = hp.npix2nside( self.npix )
        self.pixarea = hp.nside2pixarea(self.nside, degrees=True)
        self.theta, self.phi = hp.pix2ang( self.nside, np.arange(self.npix) )

        ### compute basic information statistics about map
        self.entropy     = self.base**(stats.entropy( post, base=self.base ))*self.pixarea
        self.information = self.base**(stats.information( post, base=self.base ))*self.pixarea

        ### make json representations
        jsonname = "%s_postC%s.js"%(self.label, self.tag)
        if verbose:
            print "building json represenation of postC"
            print "  "+jsonname
        self.jsPost = Json( list(stats.resample(self.postC, self.json_nside)),
                             self.output_dir,
                             self.output_url,
                             graceid    = self.graceid,
                             graceDbURL = self.graceDbURL,
                           ).saveAndUpload( jsonname )

        jsonname = "%s_cpostC%s.js"%(self.label, self.tag)
        if verbose:
            print "building json respresentation of cumulative postC"
            print "  "+jsonname

        # function resolution is somehow getting messed up here, so we resort to getattr
        self.jsCPost = Json( list(getattr(stats, '__to_cumulative')(stats.resample(self.postC, self.json_nside))),
                              self.output_dir,
                              self.output_url,
                              graceid    = self.graceid,
                              graceDbURL = self.graceDbURL,
                            ).saveAndUpload( jsonname )

    def make_mollweide(self, verbose=False):
        """
        make mollweide projections
        """
        if verbose:
            print "building mollweide projections"
        self.mollweide = dict()

        for projection, coord, post in [('astro hours mollweide', 'C', self.postC), ('mollweide', 'E', self.postE)]:

            ### generate figure
            fig, ax = mw.gen_fig_ax( self.figind, projection=projection )
            fig = Figure( fig, self.output_dir, self.output_url, graceid=self.graceid, graceDbURL=self.graceDbURL )
            self.figind += 1

            mw.heatmap( post, ax, color_map=self.color_map )
            mw.annotate( ax,
                         continents       = coord=='E',
                         continents_color = self.continents_color,
                         continents_alpha = self.continents_alpha,
                       )

            if self.transparent:
                fig.fig.patch.set_alpha(0.)
                ax.patch.set_alpha(0.)
                ax.set_alpha(0.)

            ### save just the heatmap
            figname = "%s_heatmap%s%s.%s"%(self.label, coord, self.tag, self.figtype)
            if verbose:
                print "  "+figname
            self.mollweide[coord] = fig.saveAndUpload( figname )

            ### annotate with fancy crap
            mapIND = post.argmax()
            mapY = self.theta[mapIND]
            mapX = self.phi[mapIND]
            if coord == 'C':
                mapY = 0.5*np.pi - mapY ### convert from Theta->Dec

            mw.annotate( ax,
                        projection          = projection,
                        line_of_sight       = mw.gen_line_of_sight( self.ifo_pairs, coord=coord, gps=self.gps ),
                        line_of_sight_color = self.line_of_sight_color,
                        zenith              = mw.gen_zenith( self.ifos, coord=coord, gps=self.gps ),
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
            figname = "%s_heatmap%s-annotated%s.%s"%(self.label, coord, self.tag, self.figtype)
            if verbose:
                print "  "+figname
            self.mollweide[coord+" ann"] = fig.saveAndUpload( figname )

            ### add antenna patterns as contours
            genColor = colors.getColor()
            for ind, ifo in enumerate(self.ifos):

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
                            linewidths = self.mollweide_linewidths, 
                          )
                fig.fig.text(0.01, 0.99-0.05*ind, ifo, color=color, ha='left', va='top')

            ### save heatmap + fancy crap + antenna pattern contours
            figname = "%s_heatmap%s-antennas%s.%s"%(self.label, coord, self.tag, self.figtype)
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
        obj = dict()

        for ifo1, ifo2 in self.ifo_pairs:
            ifos = "".join([ifo1, ifo2])
            if verbose:
                print "  %s - %s"%(ifo1, ifo2)

            d = dict()

            sampDt = ct.gen_sampDt( ifos, Nsamp=self.dT_Nsamp )
            maxDt = sampDt[-1]

            fig, ax = ct.genDT_fig_ax( self.figind )
            fig = Figure( fig, self.output_dir, self.output_url, graceid=self.graceid, graceDbURL=self.graceDbURL )
            self.figind += 1

            ax.set_xlim(xmin=maxDt*1e3, xmax=-maxDt*1e3) ### we work in ms here...

            if self.dT_nside:
                kde = ct.post2marg( stats.resample(self.postE, self.dT_nside), ifos, sampDt, coord='E' )
            else:
                kde = ct.post2marg( self.postE, ifos, sampDt, coord='E' )

            ### compute statistics of the marginal
            d['H'] = stats.entropy( kde, base=self.base )
            d['I'] = stats.information( kde, base=self.base )
            obj[ifos] = {'H':d['H'], 'I':d['I']}

            ### plot
            ct.plot_dT( ax, sampDt, kde, xlim_dB=self.dT_xlim_dB, color=colors.getColor().next() ) ### always use the first color!

            ### decorate
            ax.set_xlabel(r'$\Delta t_{%s}\ [\mathrm{ms}]$'%(ifos))
            ax.set_ylabel(r'$p(\Delta t_{%s}|\mathrm{data})$'%(ifos))

            ct.annotate( ax, twiny = True )

            if self.no_yticks:
                ax.set_yticklabels([])

            if self.transparent:
                fig.fig.patch.set_alpha(0.)
                ax.patch.set_alpha(0.)
                ax.set_alpha(0.)

            ### save just dT marginals
            figname = "%s_dT_%s%s.%s"%(self.label, ifos, self.tag, self.figtype)
            if verbose:
                print "    "+figname
            d['fig'] = fig.saveAndUpload( figname )

            ### annotate the plot
            ct.annotate( ax,
                         [ hp.pix2ang( self.nside, np.argmax(self.postE) ) ],
                         ifos,
                         maxDt,
                         coord   = 'E',
                         gps     = self.gps,
                         color   = self.time_delay_color,
                         alpha   = self.time_delay_alpha,
                         degrees = False,
                         twiny   = False, ### already did this
                       )

            ### save annotated dT marginals
            figname = "%s_dT_%s-annotated%s.%s"%(self.label, ifos, self.tag, self.figtype)
            if verbose:
                print "    "+figname
            d['ann fig'] = fig.saveAndUpload( figname )

            plt.close(fig.fig)
            del fig

            self.dT[ifos] = d

        ### upload json file
        jsonname = "%s_dT%s.js"%(self.label, self.tag)
        if verbose:
            print "    "+jsonname
        self.dTREF = Json( obj,
                           self.output_dir,
                           self.output_url,
                           graceid    = self.graceid,
                           graceDbURL = self.graceDbURL,
                         ).saveAndUpload( jsonname )

    def make_los(self, verbose=False):
        '''
        make line-of-sight cartesian projections and statistics
        '''
        if verbose:
            print "building line-of-sight cartesian projections"
        self.los = dict() 
        obj = dict()

        for ifo1, ifo2 in self.ifo_pairs:
            ifos = "%s%s"%(ifo1,ifo2)
            if verbose:
                print "  %s - %s"%(ifo1, ifo2)

            t, p = triangulate.line_of_sight( ifo1, ifo2, coord='E' )

            fig, ax, rproj, tproj = ct.genHist_fig_ax( self.figind )
            fig = Figure( fig, self.output_dir, self.output_url, graceid=self.graceid, graceDbURL=self.graceDbURL )
            self.figind += 1

            ### rotate
            rtheta, rphi = triangulate.rotate2pole( self.theta, self.phi, t, p )

            Nbins = max(100, int(self.npix**0.5/5))

            ### compute mutual info
            mi, Hj = triangulate.compute_mi( rtheta, rphi, Nbins, weights=self.postE )
            obj[ifos] = {'MI':mi, 'Hj':Hj}
            self.los[ifos] = {'MI':mi, 'Hj':Hj}

            ### plot
            ct.histogram2d( rtheta, 
                            rphi, 
                            ax, 
                            rproj, 
                            tproj, 
                            Nbins   = Nbins, 
                            weights = self.postE, 
                            color   = colors.getColor().next(), ### always get the first color!
                            cmap    = self.color_map 
                          )

            ### silence marginal ticks
            if self.no_yticks:
                plt.setp(rproj.get_xticklabels(), visible=False)
                plt.setp(tproj.get_yticklabels(), visible=False)

            ### save
            figname = "%s_los-%s-%s%s.%s"%(self.label, ifo1, ifo2, self.tag, self.figtype)
            if verbose:
                print "    "+figname
            self.los[ifos]['fig'] = fig.saveAndUpload( figname )

            plt.close( fig.fig )
            del fig

        ### make json
        jsonname = "%s_los%s.js"%(self.label, self.tag)
        if verbose:
            print "    "+jsonname
        self.losREF = Json( obj,
                            self.output_dir,
                            self.output_url,
                            graceid    = self.graceid,
                            graceDbURL = self.graceDbURL,
                          ).saveAndUpload( jsonname )

    def make_confidence_regions(self, verbose=False):
        '''
        compute confidence regions, statistics about them, and a plot
        we compute confidence region sizes, max(dTheta), and the number and size of modes
        '''
        if verbose:
            print "analyzing confidence regions"
        self.maxDtheta = []
        self.modes     = []
        into_modes = getattr(stats, '__into_modes') ### function resolution is getting messed, so we resort to getattr
        for cr in stats.credible_region(self.postC, self.conf):
            if len(cr):
                self.maxDtheta.append( np.arccos(stats.min_all_cos_dtheta_fast(cr, self.nside))*180/np.pi )
                self.modes.append( [self.pixarea*len(_) for _ in into_modes(self.nside, cr)] )
            else: ### no pixels in the confidence region...
                self.maxDtheta.append( 0 )
                self.modes.append( [] )

        ### write json file
        jsonname = "%s_CRStats%s.js"%(self.label, self.tag)
        if verbose:
            print "  "+jsonname
        self.jsCR = Json( {'modes':self.modes, 'maxDtheta':self.maxDtheta, 'conf':list(self.conf)},
                          self.output_dir, 
                          self.output_url, 
                          graceid    = self.graceid,
                          graceDbURL = self.graceDbURL,
                        ).saveAndUpload( jsonname )
 
        ### make confidence region figures!
        self.CR = dict()

        # size
        fig, ax = ct.genCR_fig_ax( self.figind )
        fig = Figure( fig, self.output_dir, self.output_url, graceid=self.graceid, graceDbURL=self.graceDbURL )
        self.figind += 1

        ax.semilogy( self.conf, [np.sum(_) for _ in self.modes], color=colors.getColor().next() ) ### always use the first color

        ax.set_xlim(xmin=0.0, xmax=1.0)

        ax.set_xlabel('confidence')
        ax.set_ylabel('confidence region [deg$^2$]')

        figname = "%s_CRSize%s.%s"%(self.label, self.tag, self.figtype)
        if verbose:
            print "  "+figname
        self.CR['size'] = fig.saveAndUpload( figname )

        plt.close( fig.fig )

        # max{dTheta}
        fig, ax = ct.genCR_fig_ax( self.figind )
        fig = Figure( fig, self.output_dir, self.output_url, graceid=self.graceid, graceDbURL=self.graceDbURL )
        self.figind += 1

        ax.plot( self.conf, self.maxDtheta, color=colors.getColor().next() ) ### always use the first color

        ax.set_xlim(xmin=0.0, xmax=1.0)

        ax.set_xlabel('confidence')
        ax.set_ylabel('$\max\limits_{\mathrm{CR}} \Delta\\theta$ [deg]')

        figname = "%s_CRMaxdTheta%s.%s"%(self.label, self.tag, self.figtype)
        if verbose:
            print "  "+figname
        self.CR['dTheta'] = fig.saveAndUpload( figname )

        plt.close( fig.fig )

        # num modes
        fig, ax = ct.genCR_fig_ax( self.figind )
        fig = Figure( fig, self.output_dir, self.output_url, graceid=self.graceid, graceDbURL=self.graceDbURL )
        self.figind += 1
        
        ax.plot( self.conf, [len(_) for _ in self.modes], color=colors.getColor().next() ) ### always use the first color

        ax.set_xlim(xmin=0.0, xmax=1.0)

        ax.set_xlabel('confidence')
        ax.set_ylabel('No. disjoint regions')

        figname = "%s_CRNumModes%s.%s"%(self.label, self.tag, self.figtype)
        if verbose:
            print "  "+figname
        self.CR['modes'] = fig.saveAndUpload( figname )

        plt.close( fig.fig )


    def make_antenna_patterns(self, verbose=False):
        '''
        compute antenna pattern statistics
        '''
        if verbose:
            print "computing antenna pattern statistics"
        self.ant = {}

        for ifo in self.ifos:
            Fp, Fx = detector_cache.detectors[ifo].antenna_patterns( self.theta, self.phi, 0.0 )

            mapIND = self.postE.argmax()
            self.ant[ifo] = {'map': Fp[mapIND]**2 + Fx[mapIND]**2, 'ave':np.sum(self.postE * (Fp**2 + Fx**2))}

        ### make json
        jsonname = "%s_AntStats%s.js"%(self.label, self.tag)
        if verbose:
            print "  "+jsonname
        self.jsAnt = Json( self.ant,
                           self.output_dir,
                           self.output_url,
                           graceid    = self.graceid,
                           graceDbURL = self.graceDbURL,
                         ).saveAndUpload( jsonname )

    def make_postviz(self, verbose=False ):
        """
        generate a postviz page using posterior samples
        NOT IMPLEMENTED
        """
        raise NotImplementedError('need to map posterior_samples.dat into postviz interactive html page')
        ### must write pointer into self.postviz

    def make_distanceFITS(self, verbose=False ):
        """
        generate a mapping to the distance associated with each direction in the sky and save it to a FITS file.
        should match the nside (and ordering) used in the original FITS file
        NOT IMPLEMENTED
        """
        raise NotImplementedError('need to map posteriors into "distances" and provide a FITS file with this mapping')
        ### must write pointer into self.distanceFITS = {'fits':url for fitsfile, 'C':url for mollweide proj in C, 'E':url for mollweide proj in E}
        ### this is what is expected witin self.__str__ when accessing this info

    def __str__(self):
        """
        generate html document as a string
        """
        ### set up the html doc
        doc = HTML(newlines=True)
        doc.raw_text('<!DOCTYPE html>')
        htmldoc = doc.html(lang='en')

        #----------------
        ### build header
        head = htmldoc.head()

        ### set up header for bootstrap
        head.meta(charset='utf-8')
        head.meta(contents='IE=edge')._attrs.update({'http-equiv':"X-UA-Compatible"}) ### forced into modifying _attr by "-" in attribute name
        head.meta(name="viewport", content="width=device-width, initial-scale=1")

        ### set up header for this specific template
        head.link(href="css/bootstrap.min.css", rel="stylesheet") ### need to reference full URL?
        head.link(href="starter-template.css", rel="stylesheet") ### need to reference full URL?

        ### other header information
        head.meta(name="description", content="a summary of %s"%self.fitsname)
        head.meta(name="author", content=getpass.getuser()) ### whoever ran this is the author
        
        if self.graceid:
            head.title("%s:%s"%(self.graceid, self.fitsname))
        else:
            head.title(self.fitsname)

        #----------------
        ### build body
        body = htmldoc.body()

        ### add navbar
        nav = body.nav(klass='navbar navar-inverse navebar-fixed-top') ### top level tag
        div = nav.div(klass='container')

        div1 = div.div(klass='navbar-header') ### first div, has button to take me to the top level
        link = div1.a(self.fitsname, klass='navbar-brand', href='#')

        #---
        ### add sections
        sections = body.div(klass='container')

        ### top level summary
        if hasattr(self, 'nside'): ### we must have called self.readFITS(), make_json, make_cumulative_json
            div = sections.div(id='summary', klass='container')

            row = div.div(klass='row') ### row for fitsname, nside, entropy, and information

            row.div(klass='col-md-4').a(self.fitsname, href=os.path.join(self.output_url, self.fitsname))
            row.div(klass='col-md-4').p('nside = %d'%self.nside)
            row.div(klass='col-md-4').p('H = %.3f deg2'%self.entropy)
            row.div(klass='col-md-4').p('I = %.3f deg2'%self.information)
      
            row = div.div(klass='row') ### row for probability lookup
 
            row.div(id='prob lookup', klass='col-md-3').p('RA,Dec form goes here!') ### this needs to become a field users can fill in
            row.div(id='probPerDeg2', klass='col-md-3').p('probability/deg2') ### this should reference self.jsPost
            row.div(id='cumProb', klass='col-md-3').p('cumulative probability') ### this should reference self.jsCPost
        
            if hasattr(self, 'distanceFITS'): ### must have called make_distanceFITS
                row = div.div(klass='row')

                row.div(klass='col-md-2').a('FITS file representing expected relative distances across the sky', href=self.distanceFITS['fits'])
                row.div(klass='col-md-2').img(src=self.distanceFITS['C'])                
                row.div(klass='col-md-2').img(src=self.distanceFITS['E'])              

            humanReadable = div.p('Human readable summary goes here! Should be based on the rarity of MID compared to injection sets (may be pipeline dependent...)')

        ### add mollweide plots
        if hasattr(self, 'mollweide'): ### we must have called make_mollweide
            div = sections.div(id='mollweide', klass='container')
            div.p('Mollweide projections')

            ### iterate over coordinates
            for coord in "C E".split():
                row = div.div(klass='row')

		col = row.div(klass='col-md-1')
                col.img(src=self.mollweide[coord])
                col.img(src=self.mollweide[coord+' ann'])
                col.img(src=self.mollweide[coord+' ant'])

        ### add confidence region plots
        if hasattr(self, 'CR'): ### must have called make_confidence_regions
            div = sections.div(id='confidence regions', klass='container')
            div.p('Confidence regions')

            row = div.div(klass='row')

            ### put in the figures
            col = row.div(klass='col-md-2')
            col.img(src=self.CR['size']) 
            col.img(src=self.CR['dTheta']) 
            col.img(src=self.CR['modes']) 

            ### put in the statistics

            print "WARNING: several of these should be interactive (ie: pull down) and should reference the json file, but we hack it for now"

            col = row.div(klass='col-md-2') 
            row = col.div(klass='row')
            row.div(klass='col-md-4').p('confidence')
            row.div(klass='col-md-4').p('size [deg2]')
            row.div(klass='col-md-4').p('max{dTheta} [deg]')
            row.div(klass='col-md-4').p('No. disjoint regions')

            for conf, dTheta, modes in zip(self.conf, self.maxDtheta, self.modes):
                row.div(klass='col-md-4').p('%.3f'%conf)
                row.div(klass='col-md-4').p('%.3f'%np.sum(modes))
                row.div(klass='col-md-4').p('%.3f'%dTheta)
                row.div(klass='col-md-4').p('%d'%len(modes))

        ### add antenna patterns stuff
        if hasattr(self, 'ant'): ### must have called make_antenna_patterns
            div = sections.div(id='antenna patterns', klass='container')
            div.p('Antenna Patterns')

            row = div.div(klass='row')
            row.div(klass='col-md-3').p('IFO')
            row.div(klass='col-md-3').p('F+^2 + Fx^2 @ MAP')
            row.div(klass='col-md-3').p('<F+^2 + Fx^2>')

            for ifo in self.ifos:
                row = div.div(klass='row')
                row.div(klass='col-md-3').p(ifo)
                row.div(klass='col-md-3').p('%.3f'%self.ant[ifo]['map'])
                row.div(klass='col-md-3').p('%.3f'%self.ant[ifo]['ave'])

        ### add line-of-sight sanity checks and dT marginals
        if hasattr(self, 'dT') and hasattr(self, 'los'): ### must have called make_los and make_dT
            div = sections.div(id='line of sight', klass='container')
            div.p('time delay marginals and line-of-sight frame')

            for ifo1, ifo2 in self.ifo_pairs:
                ifos = "%s%s"%(ifo1, ifo2)

                row = div.div(klass='row')

                ### first col declares ifos and gives statistics
                col = row.div(klass='col-md-3').div(klass='row').div(klass='col-md-1')
                col.p('%s - %s'%(ifo1, ifo2))
                col.p('H(dT) = %.3f'%self.dT[ifos]['H'])
                col.p('I(dT) = %.3f'%self.dT[ifos]['I'])
                col.p('MI = %.3f'%self.los[ifos]['MI'])
                col.p('Hjnt = %.3f'%self.los[ifos]['Hj'])

                ### second col contains time-delay marginals
                col = row.div(klass='col-md-3')
                with col.div(klass='row') as r:
                    r.div(klass='col-md-2').img(src=self.dT[ifos]['fig'])
                    r.div(klass='col-md-2').img(src=self.dT[ifos]['ann fig'])

                ### third col contains sanity check plot in line-of-sight frame
                row.div(klass='col-md-3').img(src=self.los[ifos]['fig'])

        ### add postviz
        if hasattr(self, 'postviz'): ### must have called make_postviz
            div = sections.div(id='postviz', klass='container')
            div.p('interactive posterior visualization')

            raise NotImplementedError('currently do not support postviz!')

        #----------------
        ### print document and return
        return str(doc)


    def write(self, verbose=False):
        '''
        writes the html document into a predictable filename
        '''
        htmlname = os.path.join( self.output_dir, "%s%s.html"%(self.label, self.tag) )
        if verbose:
            print "writing html document : "+htmlname
        file_obj = open(htmlname, "w")
        file_obj.write( str(self) )
        file_obj.close()
        return htmlname

#-------------------------------------------------

class multFITS(object):
    '''
    a class that houses data and renders html, json for comparison info about multiple FITS files

    this class also encapsulates how the html pages are structured. This should mirrored in snglFITS as much as possible: we should have "stack-posteriors" equivalents of all the plots presented in snglFITS. We may also want the user to be able to select which FITS are included in the comparison, and therefore will have to have many, many plots ready to go.
    '''

    def __init__( self, 
                  fitsnames = [],
                  output_dir='.', 
                  output_url='.', 
                  graceid=None, 
                  graceDbURL='https://gracedb.ligo.org/api/'
                ):
        self.fitsnames = fitsnames

        self.output_dir = output_dir
        self.output_url = output_url

        self.graceid    = graceid
        self.graceDbURL = graceDbURL

        ### local references for plotting
        self.figind = 0

    def __str__(self):
        raise NotImplementedError('you have a *lot* to implement before this is finished...')

    def write(self, verbose=False):
        '''
        writes the html document into a predictable filename
        '''
        htmlname = os.path.join( self.output_dir, "%s%s.html"%(self.label, self.tag) )
        if verbose:
            print "  "+htmlname
        file_obj = open(htmlname, "w")
        file_obj.write( str(self) )
        file_obj.close()
        return htmlname
