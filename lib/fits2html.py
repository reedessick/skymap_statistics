description = "a module that houses classes that write html structures for fits2html.py and friends. NOTE: there's a lot of ~repeated code between snglFITS and multFITS which could possibly be unified somehow. Something to think about..."
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

standard_tagname = ['skymapAutosummary']

#-------------------------------------------------

class Figure(object):
    '''
    a thin wrapper around figure objects that knows how to upload and save them
    '''

    def __init__(self, fig, output_dir, output_url, graceid=None, graceDbURL='https://gracedb.ligo.org/api/', upload=False):
        self.fig = fig

        self.output_dir = output_dir
        self.output_url = output_url

        self.graceid    = graceid
        self.graceDbURL = graceDbURL

        self.upload = upload

    def saveAndUpload(self, figname, message='', tagname=[]):
        filename = os.path.join(self.output_dir, figname)
        self.fig.savefig( filename )
        if self.upload and (self.graceid!=None):
            gdb = GraceDb(self.graceDbURL)
            httpResponse = gdb.writeLog( self.graceid, message=message, filename=filename, tagname=standard_tagname+tagname )
            ### may want to check httpResponse for errors...

        return os.path.join(self.output_url, figname) 

class Json(object):
    '''
    a thin wrapper around a json object that knows how to upload and save them
    '''

    def __init__(self, obj, output_dir, output_url, graceid=None, graceDbURL='https://gracedb.ligo.org/api/', upload=False):
        self.obj = obj

        self.output_dir = output_dir
        self.output_url = output_url

        self.graceid    = graceid
        self.graceDbURL = graceDbURL

        self.upload = upload

    def saveAndUpload(self, filename, message='', tagname=[]):
        fileName = os.path.join(self.output_dir, filename)
        file_obj = open( fileName, "w" )
        json.dump( self.obj, file_obj )
        file_obj.close()
        if self.upload and (self.graceid!=None):
            gdb = GraceDb(self.graceDbURL)
            httpResponse = gdb.writeLog( self.graceid, message=message, filename=fileName, tagname=standard_tagname+tagname )
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
                  upload     = False,
                  ### options for json reference files
                  json_nside = 128,
                  ### general options about annotation and which plots to build
                  ifos = [],
                  ### general options about colors, shading, and labeling
                  color_map    = "OrRd",
                  transparent  = False,
                  no_margticks = False,
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
        self.label    = os.path.basename(fitsname)
        if self.label.endswith('.gz'):
            self.label = self.label[:-3] ### get rid of .gz
        self.label = self.label[:-5] ### get rid of .fits

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
        self.upload     = upload

        ### general color schemes
        self.color_map    = color_map
        self.transparent  = transparent
        self.no_margticks = no_margticks

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
            print "reading : %s -> %s"%(self.fitsname, self.label)
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
        self.npix    = len(post)
        self.nside   = hp.npix2nside( self.npix )
        self.pixarea = hp.nside2pixarea(self.nside, degrees=True)
        self.theta, self.phi = hp.pix2ang( self.nside, np.arange(self.npix) )

        ### compute basic information statistics about map
        self.entropy     = stats.entropy( post, base=self.base )
        self.information = stats.information( post, base=self.base )

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
                             upload     = self.upload,
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
                              upload     = self.upload,
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
            fig = Figure( fig, self.output_dir, self.output_url, graceid=self.graceid, graceDbURL=self.graceDbURL, upload=self.upload )
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
            for ind, ifo in enumerate(self.ifos):

                Fp, Fx = detector_cache.detectors[ifo].antenna_patterns( self.theta, self.phi, 0.0 )
                ant = Fp**2 + Fx**2
                ant /= np.sum(ant)
                if coord == 'C':
                    ant = triangulate.rotateMapE2C( ant, self.gps )

                color = colors.getIFOColor( ifo )
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

            ### build contour plot
            fig, ax = mw.gen_fig_ax( self.figind, projection=projection )
            fig = Figure( fig, self.output_dir, self.output_url, graceid=self.graceid, graceDbURL=self.graceDbURL, upload=self.upload )
            self.figind += 1

            mw.contour( post, 
                        ax, 
                        colors     = colors.getColor().next(), ### always use the first color!
                        levels     = self.mollweide_levels,
                        alpha      = self.mollweide_alpha,
                        linewidths = self.mollweide_linewidths,
                      )
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
            figname = "%s_contour%s%s.%s"%(self.label, coord, self.tag, self.figtype)
            if verbose:
                print "  "+figname
            self.mollweide[coord+" cnt"] = fig.saveAndUpload( figname )

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
            fig = Figure( fig, self.output_dir, self.output_url, graceid=self.graceid, graceDbURL=self.graceDbURL, upload=self.upload )
            self.figind += 1

            ax.set_xlim(xmin=maxDt*1e3, xmax=-maxDt*1e3) ### we work in ms here...

            if self.dT_nside:
                kde = ct.post2marg( stats.resample(self.postE, self.dT_nside), ifos, sampDt, coord='E' )
            else:
                kde = ct.post2marg( self.postE, ifos, sampDt, coord='E' )

            ### compute statistics of the marginal
            d['H'] = stats.entropy( kde, base=self.base )
            d['I'] = stats.information( kde, base=self.base )
            d['thetaMAP'] = np.arccos( sampDt[kde.argmax()]/maxDt )
            obj[ifos] = {'H':d['H'], 'I':d['I'], 'thetaMAP':d['thetaMAP']}

            ### plot
            ct.plot_dT( ax, sampDt, kde, xlim_dB=self.dT_xlim_dB, color=colors.getColor().next() ) ### always use the first color!

            ### decorate
            ax.set_xlabel(r'$\Delta t_{%s}\ [\mathrm{ms}]$'%(ifos))
            ax.set_ylabel(r'$p(\Delta t_{%s}|\mathrm{data})$'%(ifos))

            ct.annotate( ax, IFOs=ifos, twiny=True )

            if self.no_margticks:
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
                           upload     = self.upload,
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
            fig = Figure( fig, self.output_dir, self.output_url, graceid=self.graceid, graceDbURL=self.graceDbURL, upload=self.upload )
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
            if self.no_margticks:
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
                            upload     = self.upload,
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
                          upload     = self.upload,
                        ).saveAndUpload( jsonname )
 
        ### make confidence region figures!
        self.CR = dict()

        # size
        fig, ax = ct.genCR_fig_ax( self.figind )
        fig = Figure( fig, self.output_dir, self.output_url, graceid=self.graceid, graceDbURL=self.graceDbURL, upload=self.upload )
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
        fig = Figure( fig, self.output_dir, self.output_url, graceid=self.graceid, graceDbURL=self.graceDbURL, upload=self.upload )
        self.figind += 1

        ax.plot( self.conf, self.maxDtheta, color=colors.getColor().next() ) ### always use the first color

        ax.set_xlim(xmin=0.0, xmax=1.0)

        ax.set_xlabel('confidence')
        ax.set_ylabel('$\max\limits_{\mathrm{CR}}\left\{ \Delta\\theta \\right\}$ [deg]')

        figname = "%s_CRMaxdTheta%s.%s"%(self.label, self.tag, self.figtype)
        if verbose:
            print "  "+figname
        self.CR['dTheta'] = fig.saveAndUpload( figname )

        plt.close( fig.fig )

        # num modes
        fig, ax = ct.genCR_fig_ax( self.figind )
        fig = Figure( fig, self.output_dir, self.output_url, graceid=self.graceid, graceDbURL=self.graceDbURL, upload=self.upload )
        self.figind += 1
        
        ax.plot( self.conf, [len(_) for _ in self.modes], color=colors.getColor().next() ) ### always use the first color

        ax.set_xlim(xmin=0.0, xmax=1.0)
        ax.set_ylim(ymin=0.0)
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
                           upload     = self.upload,
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

    def __str__(self):
        """
        generate html document as a string
        """
        fitsbasename = os.path.basename(self.fitsname)

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
        head.link( href        = "https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css", 
                   rel         = "stylesheet",
                   integrity   = "sha384-BVYiiSIFeK1dGmJRAkycuHAHRg32OmUcww7on3RYdg4Va+PmSTsz/K68vbdEjh4u", 
                   crossorigin = "anonymous",
                 )

        ### other header information
        head.meta(name="description", content="a summary of %s"%fitsbasename)
        head.meta(name="author", content=getpass.getuser()) ### whoever ran this is the author
        
        if self.graceid:
            head.title("%s:%s"%(self.graceid, fitsbasename))
        else:
            head.title(fitsbasename)

        #----------------
        ### build body
        body = htmldoc.body()

        ### add navbar
        div1 = body.nav(klass='navbar navar-inverse navebar-fixed-top').div(klass='container').div(klass='navbar-header') ### first div, has links to sections
        div1.a(self.label, klass='navbar-brand')

        div1 = body.nav(klass='navbar navar-inverse navebar-fixed-top').div(klass='container').div(klass='navbar-header') ### second div, contains links to sections

        #---
        ### add sections
        sections = body.div(klass='container')

        ### top level summary
        if hasattr(self, 'nside'): ### we must have called self.readFITS(), make_json, make_cumulative_json
            sections.hr
            div = sections.div(id='summary', klass='container')
            div.h1('Summary', id='summary')
            div1.a('Summary', klass='navbar-brand', href='#summary') ### add to top-level navigation bar

            row = div.div(klass='row') ### row for fitsname, nside, entropy, and information

            row.div(klass='col-md-3').a(fitsbasename, href=os.path.join(self.output_url, fitsbasename))
            row.div(klass='col-md-2').p('nside = %d'%self.nside)
            row.div(klass='col-md-3').p().raw_text('H = %.3f (%.3f deg<sup>2</sup>)'%(self.entropy, self.base**(self.entropy)*self.pixarea))
            row.div(klass='col-md-2').p('I = %.3f'%self.information)
     
            ### give maximum a posteriori positions
            dec, ra = np.array(hp.pix2ang(self.nside, self.postC.argmax()))*180/np.pi
            dec = 90 - dec
            div.div(klass='row').div(klass='col-md-5').p().raw_text('maximum a posteriori (&alpha;, &delta;) = (%.2f&deg;, %.2f&deg;)'%(ra, dec))

            lat, lon = np.array(hp.pix2ang(self.nside, self.postE.argmax()))*180/np.pi
            lat = 90 - lat
            div.div(klass='row').div(klass='col-md-5').p().raw_text('maximum a posteriori (lon, lat) = (%.2f&deg;, %.2f&deg;)'%(lon, lat))
 
            ### row for probability lookup
            row = div.div(klass='row')
 
            row.div(id='prob lookup', klass='col-md-3').p().b().raw_text('&alpha;, &delta; form goes here!') ### this needs to become a field users can fill in
            row.div(id='probPerDeg2', klass='col-md-2').p().raw_text('probability/deg<sup>2</sup>') ### this should reference self.jsPost
            row.div(id='cumProb', klass='col-md-2').p('cumulative probability') ### this should reference self.jsCPost
      
            if hasattr(self, 'distanceFITS'): ### must have called make_distanceFITS
                row = div.div(klass='row')

                row.div(klass='col-md-2').a('FITS file representing expected relative distances across the sky', href=self.distanceFITS['fits'])
                row.div(klass='col-md-2').img(src=self.distanceFITS['C'])                
                row.div(klass='col-md-2').img(src=self.distanceFITS['E'])              

            ### add human readable summary
            humanReadable = div.p()
#            humanReadable.b('Human readable summary goes here!') 

            if hasattr(self, 'los'): ### statements about the relative rarety of the observed mutual information distance
                humanReadable.br
                mid = np.array([self.los["%s%s"%ifos]['MI']/self.los["%s%s"%ifos]['Hj'] for ifos in self.ifo_pairs])
                if np.any(mid > 0.1):
                    humanReadable += 'Mutual Information Distance (MID) is <b>large</b> (> 0.1) for at least one IFO pair!'
                else:
                    humanReadable += 'Mutual Information Distance (MID) is <b>small</b> (< 0.1) for all IFO pairs.'

        ### add mollweide plots
        if hasattr(self, 'mollweide'): ### we must have called make_mollweide
            sections.hr
            div = sections.div(id='mollweide', klass='container')
            div.h1('Mollweide Projections', id='mollweide')
            div1.a('Mollweide Projections', klass='navbar-brand', href='#mollweide') ### add to top-level navigation bar

            ### iterate over coordinates
            width = '575' ### FIXME: hard coded figure width isn't great...
            for coord, label in [("C","Equatorial"), ("E", "Geographic")]:
                div.h2(label+' coordinates', id="mollweide"+label)

                row = div.div(klass='row')
                row.div(klass='col-md-6').img(src=self.mollweide[coord], width=width)
                row.div(klass='col-md-6').img(src=self.mollweide[coord+' cnt'], width=width)

                row = div.div(klass='row')
                row.div(klass='col-md-6').img(src=self.mollweide[coord+' ann'], width=width)
                row.div(klass='col-md-6').img(src=self.mollweide[coord+' ant'], width=width)

        ### add confidence region plots
        if hasattr(self, 'CR'): ### must have called make_confidence_regions
            sections.hr
            div = sections.div(id='confidence regions', klass='container')
            div.h1('Confidence Regions', id='confidence')
            div1.a('Confidence Regions', klass='navbar-brand', href='#confidence') ### add to top-level navigation bar

            row = div.div(klass='row')

            ### put in the figures
            width = '550' ### FIXME: hard coding width isn't great...
            col = row.div(klass='col-md-6')
            col.img(src=self.CR['size'], width=width) 
            col.img(src=self.CR['dTheta'], width=width) 
            col.img(src=self.CR['modes'], width=width) 

            ### put in the statistics

            print "\nWARNING: several of these should be interactive (ie: pull down) and should reference the json file, but we hack it for now\n"

            col = row.div(klass='col-md-6') 
            row = col.div(klass='row')
            row.div(klass='col-md-2').p('confidence', align='center')
            row.div(klass='col-md-2').p(align='right').raw_text('size [deg<sup>2</sup>]')
            row.div(klass='col-md-3').p(align='right').raw_text('max{&Delta;&theta;} [&deg;]')
            row.div(klass='col-md-4').p('No. disjoint regions', align='center')

            for conf, dTheta, modes in zip(self.conf, self.maxDtheta, self.modes):
                row = col.div(klass='row')
                row.div(klass='col-md-2').p('%.1f%s'%(100*conf, "%"), align='right')
                row.div(klass='col-md-2').p('%.3f'%np.sum(modes), align='right')
                row.div(klass='col-md-3').p('%.3f'%dTheta, align='right')
                row.div(klass='col-md-3').p('%d'%len(modes), align='right')

        ### add antenna patterns stuff
        if hasattr(self, 'ant'): ### must have called make_antenna_patterns
            sections.hr
            div = sections.div(id='antenna patterns', klass='container')
            div.h1('Antenna Patterns', id='antenna')
            div1.a('Antenna Patterns', klass='navbar-brand', href='#antenna') ### add to top-level navigation bar

            row = div.div(klass='row')
            row.div(klass='col-md-1').p('IFO', align='center')
            row.div(klass='col-md-2').p(align='center').raw_text('(F<sub>+</sub><sup>2</sup> + F<sub>x</sub><sup>2</sup>)<sup>1/2</sup> @ MAP')
            row.div(klass='col-md-2').p(align='center').raw_text('mean{F<sub>+</sub><sup>2</sup> + F<sub>x</sub><sup>2</sup>}<sup>1/2</sup>')

            for ifo in self.ifos:
                row = div.div(klass='row')
                row.div(klass='col-md-1').p(ifo, align='center')
                row.div(klass='col-md-2').p('%.3f'%(self.ant[ifo]['map']**0.5), align='center')
                row.div(klass='col-md-2').p('%.3f'%(self.ant[ifo]['ave']**0.5), align='center')

        ### add line-of-sight sanity checks and dT marginals
        if hasattr(self, 'dT') and hasattr(self, 'los'): ### must have called make_los and make_dT
            sections.hr
            div = sections.div(id='line of sight', klass='container')
            div.h1('Time Delay Marginals and Line-of-Sight Frame', id='timeDelay')
            div1.a('Time Delay', klass='navbar-brand', href='#timeDelay') ### add to top-level navigation bar

            for ifo1, ifo2 in self.ifo_pairs:
                div.h2('%s - %s'%(ifo1, ifo2))
                row = div.div(klass='row')

                ifos = "%s%s"%(ifo1, ifo2)
                ### first col declares ifos and gives statistics
                col = row.div(klass='col-md-2')#.div(klass='row').div(klass='col-md-3')
                col.p('H(dT) = %.3f'%self.dT[ifos]['H'], align='center')
                col.p('I(dT) = %.3f'%self.dT[ifos]['I'], align='center')
                col.p(align='center').raw_text('&theta;<sub>MAP</sub> = %.2f&deg;'%(self.dT[ifos]['thetaMAP']*180/np.pi))
                col.p('MI = %.3f'%self.los[ifos]['MI'], align='center')
                col.p(align='center').raw_text('H<sub>jnt</sub> = %.3f'%self.los[ifos]['Hj'])
                col.p('MID = %.5f'%(self.los[ifos]['MI']/self.los[ifos]['Hj']), align='center')

                ### second col contains time-delay marginals
                col = row.div(klass='col-md-8')
                width = '700' ### FIXME: hard coding width isn't great...
#               col.img(src=self.dT[ifos]['fig'], width=width)
                col.img(src=self.dT[ifos]['ann fig'], width=width)
                width = '900' ### FIXME: hard coding width isn't great...
                col.img(src=self.los[ifos]['fig'], width=width)

        ### add postviz
        if hasattr(self, 'postviz'): ### must have called make_postviz
            sections.hr
            div = sections.div(id='postviz', klass='container')
            div.h1('Interactive Posterior Visualization', id='postviz')
            div1.a('Postviz', klass='navbar-brand', href='#postviz') ### add to top-level navigation bar

            raise NotImplementedError('currently do not support postviz!')

        #----------------
        ### print document and return
        return str(doc)

#-------------------------------------------------

class multFITS(object):
    '''
    a class that houses data and renders html, json for comparison info about multiple FITS files

    this class also encapsulates how the html pages are structured. This should mirrored in snglFITS as much as possible: we should have "stack-posteriors" equivalents of all the plots presented in snglFITS. We may also want the user to be able to select which FITS are included in the comparison, and therefore will have to have many, many plots ready to go.
    '''

    def __init__( self, 
                  fitsnames,
                  ### general options about output
                  output_dir = '.', 
                  output_url = './', ### ignored if graceid is supplied, otherwise used to build URLs in html document
                  tag        = '',
                  figtype    = "png",
                  dpi        = 500,
                  graceid    = None, ### if supplied, upload files and reference them in the html document
                  graceDbURL = 'https://gracedb.ligo.org/api/',
                  upload     = False,
                  ### general options about annotation and which plots to build
                  ifos = [],
                  ### general options about colors, shading, and labeling
                  color_map    = "OrRd",
                  transparent  = False,
                  no_margticks = False,
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
                  conf = np.linspace(0, 1.0, 51),
                  area = np.logspace(1, 4, 51),
                ):

        ### general things about FITS file
        self.fitsnames = sorted(fitsnames)
        self.fits_pairs = []
        for ind, fits1 in enumerate(self.fitsnames):
            for fits2 in self.fitsnames[ind+1:]:
                self.fits_pairs.append( (fits1, fits2) )

#        self.labels = dict((fitsname, os.path.basename(fitsname).split('.')[0]) for fitsname in fitsnames)
        self.labels = {}
        for fitsname in fitsnames:
            label = os.path.basename(fitsname)
            if label.endswith('.gz'): ### get rid of .gz
                label = label[:-3]
            label = label[:-5] ### get rid of .fits
            self.labels[fitsname] = label

        self.texlabels = dict((key, val.replace('_','\_')) for key, val in self.labels.items())

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
            self.label = graceid
        else:
            self.output_url = output_url
            self.label = 'multFITS'

        self.tag   = tag

        self.figtype = figtype
        self.dpi     = dpi

        self.graceid    = graceid
        self.graceDbURL = graceDbURL
        self.upload     = upload

        ### general color schemes
        self.color_map    = color_map
        self.transparent  = transparent
        self.no_margticks = no_margticks

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
        self.area = area

        ### local references for plotting
        self.figind = 0

    def readFITS(self, verbose=False):
        '''
        reads in the FITS file and sets up local copies
        '''
        self.fitsdata = {}
        for fitsname in self.fitsnames:
            data = {}
            if verbose:
                print "reading : %s -> %s"%(fitsname, self.labels[fitsname])
            ### load in map
            post, header = hp.read_map( fitsname, h=True, verbose=verbose )
            header = dict(header)

            ### ensure we are in RING ordering
            if header['ORDERING']=='NEST':
                post = hp.reorder( post, 'NEST', r2n=1 )

            ### extract gps time
            data['gps'] = tconvert(header['DATE-OBS'])

            ### set up references to maps in C and E coordinates
            if verbose:
                print "  setting up local copies in C and E coordinates"
            coord = header['COORDSYS']
            if coord == 'C':
                data['C'] = postC = post[:]
                data['E'] = triangulate.rotateMapC2E( data['C'], data['gps'] )
            elif coord == 'E':
                data['E'] = post[:]
                data['C'] = triangulate.rotateMapE2C( data['E'], data['gps'] )
            else:
                raise ValueError('COORDSYS=%s not understood!'%coord)

            ### set up meta-data about the map
            if verbose:
                print "  setting up local references to angles"
            data['npix']    = len(post)
            data['nside']   = hp.npix2nside( data['npix'] )
            data['pixarea'] = hp.nside2pixarea(data['nside'], degrees=True)

            ### add to local structure
            self.fitsdata[fitsname] = data

    def make_mollweide(self, verbose=False):
        """
        make mollweide projections
        """
        if verbose:
            print "building mollweide projections"
        self.mollweide = dict()

        for projection, coord in [('astro hours mollweide', 'C'), ('mollweide', 'E')]:

            ### generate figure
            fig, ax = mw.gen_fig_ax( self.figind, projection=projection )
            fig = Figure( fig, self.output_dir, self.output_url, graceid=self.graceid, graceDbURL=self.graceDbURL, upload=self.upload )
            self.figind += 1

            getColor = colors.getColor()
            for ind, fitsname in enumerate(self.fitsnames): ### iterate through FITS files
                if verbose:
                    print "    "+self.labels[fitsname]
                color = getColor.next()
                mw.contour( self.fitsdata[fitsname][coord],
                            ax,
                            colors     = color,
                            levels     = self.mollweide_levels,
                            alpha      = self.mollweide_alpha,
                            linewidths = self.mollweide_linewidths,
                          )
                fig.fig.text(0.01, 0.99-0.05*(ind), self.texlabels[fitsname], color=color, ha='left', va='top')
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
            figname = "%s_contour%s%s.%s"%(self.label, coord, self.tag, self.figtype)
            if verbose:
                print "  "+figname
            self.mollweide[coord] = fig.saveAndUpload( figname )

            ### annotate with fancy crap
            getColor = colors.getColor()
            for ind, fitsname in enumerate(self.fitsnames):
                if verbose:
                    print "    "+self.labels[fitsname]
                mapIND = self.fitsdata[fitsname][coord].argmax()
                theta, phi = hp.pix2ang( self.fitsdata[fitsname]['nside'], np.arange(self.fitsdata[fitsname]['npix']) )
                mapY = theta[mapIND]
                mapX = phi[mapIND]
                if coord == 'C':
                    mapY = 0.5*np.pi - mapY ### convert from Theta->Dec

                gps = self.fitsdata[fitsname]['gps']
                color = getColor.next()

                mw.annotate( ax,
                            projection          = projection,
                            line_of_sight       = mw.gen_line_of_sight( self.ifo_pairs, coord=coord, gps=gps ),
                            line_of_sight_color = color,
                            zenith              = mw.gen_zenith( self.ifos, coord=coord, gps=gps ),
                            zenith_color        = color,
                            time_delay          = mw.gen_time_delay( [(mapY, mapX)], self.ifo_pairs, coord=coord, gps=gps, degrees=False ),
                            time_delay_color    = color,
                            time_delay_alpha    = self.time_delay_alpha,
                            marker_Dec_RA       = mw.gen_marker_Dec_RA( [(mapY, mapX)], coord=coord, gps=gps, degrees=False ),
                            marker              = self.marker,
                            marker_color        = color,
                            marker_size         = self.marker_size,
                            marker_edgewidth    = self.marker_edgewidth,
                            marker_alpha        = self.marker_alpha,
                          )

            ### save heatmap + fancy crap
            figname = "%s_contour%s-annotated%s.%s"%(self.label, coord, self.tag, self.figtype)
            if verbose:
                print "  "+figname
            self.mollweide[coord+" ann"] = fig.saveAndUpload( figname )

    def make_dT(self, verbose=False):
        '''
        make time-delay marginal plots and statistics
        '''
        if verbose:
            print "building time-delay marginals"
        self.dT = dict()
        obj = {}

        for ifo1, ifo2 in self.ifo_pairs:
            ifos = "".join([ifo1, ifo2])
            if verbose:
                print "  %s - %s"%(ifo1, ifo2)

            d = dict()
            dd = dict()

            sampDt = ct.gen_sampDt( ifos, Nsamp=self.dT_Nsamp )
            maxDt = sampDt[-1]

            fig, ax = ct.genDT_fig_ax( self.figind )
            fig = Figure( fig, self.output_dir, self.output_url, graceid=self.graceid, graceDbURL=self.graceDbURL, upload=self.upload )
            self.figind += 1

            ax.set_xlim(xmin=maxDt*1e3, xmax=-maxDt*1e3) ### we work in ms here...

            getColor = colors.getColor() ### restart for each posterior
            margs = {}
            for ind, fitsname in enumerate(self.fitsnames):
                if verbose:
                    print "      "+self.labels[fitsname]

                print "\nWARNING: may want to read in the json file with marginal distribution instead of recomputing it...\n"

                if self.dT_nside:
                    kde = ct.post2marg( stats.resample(self.fitsdata[fitsname]['E'], self.dT_nside), ifos, sampDt, coord='E' )
                else:
                    kde = ct.post2marg( self.fitsdata[fitsname]['E'], ifos, sampDt, coord='E' )

                margs[fitsname] = kde
                dd[fitsname] = {'thetaMAP' : np.arccos( sampDt[kde.argmax()]/maxDt )}

                ### plot
                color = getColor.next()
                ct.plot_dT( ax, sampDt, kde, xlim_dB=self.dT_xlim_dB, color=color )
                fig.fig.text(0.10+0.02, 0.93-0.05*ind, self.texlabels[fitsname], color=color, ha='left', va='top')

            ### compute similarity measures
            for fits1, fits2 in self.fits_pairs:
                dd["%s|%s"%(fits1, fits2)] = {'fidelity'  : stats.fidelity( margs[fits1], margs[fits2] ) }
            d.update( dd )
            obj[ifos] = dd

            ### decorate
            ax.set_xlabel(r'$\Delta t_{%s}\ [\mathrm{ms}]$'%(ifos))
            ax.set_ylabel(r'$p(\Delta t_{%s}|\mathrm{data})$'%(ifos))

            ct.annotate( ax, IFOs=ifos, twiny=True )

            if self.no_margticks:
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
            getColor = colors.getColor()
            for fitsname in self.fitsnames:
                if verbose:
                    print "      "+self.labels[fitsname]
                ct.annotate( ax,
                             [ hp.pix2ang( self.fitsdata[fitsname]['nside'], np.argmax(self.fitsdata[fitsname]['E']) ) ],
                             ifos,
                             coord   = 'E',
                             gps     = self.fitsdata[fitsname]['gps'],
                             color   = getColor.next(),
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
                           upload     = self.upload,
                         ).saveAndUpload( jsonname )


    def make_los(self, verbose=False):
        '''
        make line-of-sight cartesian projections and statistics
        '''
        if verbose:
            print "building line-of-sight cartesian projections"
        self.los = dict()

        for ifo1, ifo2 in self.ifo_pairs:
            ifos = "%s%s"%(ifo1,ifo2)
            if verbose:
                print "  %s - %s"%(ifo1, ifo2)

            t, p = triangulate.line_of_sight( ifo1, ifo2, coord='E' )

            fig, ax, rproj, tproj = ct.genHist_fig_ax( self.figind )
            fig = Figure( fig, self.output_dir, self.output_url, graceid=self.graceid, graceDbURL=self.graceDbURL, upload=self.upload )
            self.figind += 1

            ### rotate
            getColor = colors.getColor()
            for ind, fitsname in enumerate(self.fitsnames):
                if verbose:
                    print "      "+self.labels[fitsname]
                theta, phi = hp.pix2ang( self.fitsdata[fitsname]['nside'], np.arange(self.fitsdata[fitsname]['npix']) )
                rtheta, rphi = triangulate.rotate2pole( theta, phi, t, p )

                Nbins = max(100, int(self.fitsdata[fitsname]['npix']**0.5/5))

                ### plot
                color = getColor.next()
                ct.histogram2d( rtheta,
                                rphi,
                                ax,
                                rproj,
                                tproj,
                                Nbins   = Nbins,
                                weights = self.fitsdata[fitsname]['E'],
                                color   = color,
                                contour = True,
                                alpha   = self.mollweide_alpha, ### may want to allow a separate option, but this should do for now... 
                              )
                fig.fig.text(0.95, 0.95-0.05*ind, self.texlabels[fitsname], color=color, ha='right', va='top')

            ### silence marginal ticks
            if self.no_margticks:
                plt.setp(rproj.get_xticklabels(), visible=False)
                plt.setp(tproj.get_yticklabels(), visible=False)

            ### save
            figname = "%s_los-%s-%s%s.%s"%(self.label, ifo1, ifo2, self.tag, self.figtype)
            if verbose:
                print "    "+figname
            self.los[ifos] = fig.saveAndUpload( figname )

            plt.close( fig.fig )

    def make_confidence_regions(self, verbose=False):
        '''
        compute confidence regions, statistics about them, and a plot
        we compute confidence region sizes, max(dTheta), and the number and size of modes
        NOTE: we compute geometric overlaps of confidence regions in make_comparison because that standardized the location in which we upsample maps for comparisons
        '''
        if verbose:
            print "analyzing confidence regions"

        ### make confidence region figures!
        self.CR = dict()

        # size
        fig, sizax = ct.genCR_fig_ax( self.figind )
        sizfig = Figure( fig, self.output_dir, self.output_url, graceid=self.graceid, graceDbURL=self.graceDbURL, upload=self.upload )
        self.figind += 1

        # max{dTheta}
        fig, mdtax = ct.genCR_fig_ax( self.figind )
        mdtfig = Figure( fig, self.output_dir, self.output_url, graceid=self.graceid, graceDbURL=self.graceDbURL, upload=self.upload )
        self.figind += 1

        # num modes
        fig, numax = ct.genCR_fig_ax( self.figind )
        numfig = Figure( fig, self.output_dir, self.output_url, graceid=self.graceid, graceDbURL=self.graceDbURL, upload=self.upload )
        self.figind += 1

        ### compute confidence regions and add to figure
        into_modes = getattr(stats, '__into_modes') ### function resolution is getting messed, so we resort to getattr

        getColor = colors.getColor()
        for ind, fitsname in enumerate(self.fitsnames):
            if verbose:
                print "    "+self.labels[fitsname]

            print "\nWARNING: may want to read in the json file with confidence regions instead of recomputing them...\n"

            maxDtheta = []
            modes     = []
            nside   = self.fitsdata[fitsname]['nside']
            pixarea = self.fitsdata[fitsname]['pixarea']
            for cr in stats.credible_region(self.fitsdata[fitsname]['C'], self.conf):
                if len(cr):
                    maxDtheta.append( np.arccos(stats.min_all_cos_dtheta_fast(cr, nside))*180/np.pi )
                    modes.append( [pixarea*len(_) for _ in into_modes(nside, cr)] )
                else: ### no pixels in the confidence region...
                    maxDtheta.append( 0 )
                    modes.append( [] )

            color = getColor.next()

            sizax.semilogy( self.conf, [np.sum(_) for _ in modes], color=color )
            sizfig.fig.text(0.10+0.02, 0.93-0.05*ind, self.texlabels[fitsname], color=color, ha='left', va='top')

            mdtax.plot( self.conf, maxDtheta, color=color )
            mdtfig.fig.text(0.10+0.02, 0.93-0.05*ind, self.texlabels[fitsname], color=color, ha='left', va='top')

            numax.plot( self.conf, [len(_) for _ in modes], color=color )
            numfig.fig.text(0.10+0.02, 0.93-0.05*ind, self.texlabels[fitsname], color=color, ha='left', va='top')

        # size
        sizax.set_xlim(xmin=0.0, xmax=1.0)

        sizax.set_xlabel('confidence')
        sizax.set_ylabel('confidence region [deg$^2$]')

        figname = "%s_CRSize%s.%s"%(self.label, self.tag, self.figtype)
        if verbose:
            print "  "+figname
        self.CR['size'] = sizfig.saveAndUpload( figname )

        plt.close( sizfig.fig )

        # max{dTheta}
        mdtax.set_xlim(xmin=0.0, xmax=1.0)

        mdtax.set_xlabel('confidence')
        mdtax.set_ylabel('$\max\limits_{\mathrm{CR}}\left\{ \Delta\\theta \\right\}$ [deg]')

        figname = "%s_CRMaxdTheta%s.%s"%(self.label, self.tag, self.figtype)
        if verbose:
            print "  "+figname
        self.CR['dTheta'] = mdtfig.saveAndUpload( figname )

        plt.close( mdtfig.fig )

        # num modes
        numax.set_xlim(xmin=0.0, xmax=1.0)
        numax.set_ylim(ymin=0.0)
        numax.set_xlabel('confidence')
        numax.set_ylabel('No. disjoint regions')

        figname = "%s_CRNumModes%s.%s"%(self.label, self.tag, self.figtype)
        if verbose:
            print "  "+figname
        self.CR['modes'] = numfig.saveAndUpload( figname )

        plt.close( numfig.fig )

    def make_comparison(self, verbose=False):
        '''
        generates comparison statistics using the full maps. 
        NOTE: we re-compute confidence regions here (possibly expensive!) because we need to upsample maps if their nsides differ.
              the only comparison statistics that are not copmuted here involve the time-delay marginals, which are computed in make_dT
        '''
        if verbose:
            print "computing comparison statistics"
        self.comp = dict()

        # conf region intersection
        fig, cri_ax = ct.genCR_fig_ax( self.figind )
        cri_fig = Figure( fig, self.output_dir, self.output_url, graceid=self.graceid, graceDbURL=self.graceDbURL, upload=self.upload )
        self.figind += 1

        # conf region union
        fig, cru_ax = ct.genCR_fig_ax( self.figind )
        cru_fig = Figure( fig, self.output_dir, self.output_url, graceid=self.graceid, graceDbURL=self.graceDbURL, upload=self.upload )
        self.figind += 1

        # conf region ratio
        fig, crr_ax = ct.genCR_fig_ax( self.figind )
        crr_fig = Figure( fig, self.output_dir, self.output_url, graceid=self.graceid, graceDbURL=self.graceDbURL, upload=self.upload )
        self.figind += 1

        # conf region contained
        fig, crc_ax = ct.genCR_fig_ax( self.figind )
        crc_tx = crc_ax.twinx()
        crc_fig = Figure( fig, self.output_dir, self.output_url, graceid=self.graceid, graceDbURL=self.graceDbURL, upload=self.upload )
        self.figind += 1

        # area intersection
        fig, ari_ax = ct.genCR_fig_ax( self.figind )
        ari_fig = Figure( fig, self.output_dir, self.output_url, graceid=self.graceid, graceDbURL=self.graceDbURL, upload=self.upload )
        self.figind += 1

        # area union
        fig, aru_ax = ct.genCR_fig_ax( self.figind )
        aru_fig = Figure( fig, self.output_dir, self.output_url, graceid=self.graceid, graceDbURL=self.graceDbURL, upload=self.upload )
        self.figind += 1

        # area ratio
        fig, arr_ax = ct.genCR_fig_ax( self.figind )
        arr_fig = Figure( fig, self.output_dir, self.output_url, graceid=self.graceid, graceDbURL=self.graceDbURL, upload=self.upload )
        self.figind += 1

        # area contained
        fig, arc_ax = ct.genCR_fig_ax( self.figind )
        arc_tx = arc_ax.twinx()
        arc_fig = Figure( fig, self.output_dir, self.output_url, graceid=self.graceid, graceDbURL=self.graceDbURL, upload=self.upload )
        self.figind += 1

        ### set up labeling for contained plots
        crc_fig.fig.text(0.50, 0.93, 'A - B', color='k', ha='center', va='top')
        arc_fig.fig.text(0.50, 0.93, 'A - B', color='k', ha='center', va='top')

        ### iterate through pairs
        getColor = colors.getColor()
        for ind, (fits1, fits2) in enumerate(self.fits_pairs):

            ### upsample as needed
            nside = max( self.fitsdata[fits1]['nside'], self.fitsdata[fits2]['nside'] )
            pixarea = hp.nside2pixarea( nside, degrees=True ) ### figure out indecies for area overlaps

            post1 = stats.resample( self.fitsdata[fits1]['C'][:], nside )
            post2 = stats.resample( self.fitsdata[fits2]['C'][:], nside )

            post1argsort = post1.argsort()[::-1] ### reverse so big stuff comes first
            post2argsort = post2.argsort()[::-1]

            t1, p1 = hp.pix2ang( nside, post1.argmax() )
            t2, p2 = hp.pix2ang( nside, post2.argmax() )

            ### things involving confidence regions
            conf_I = []
            conf_U = []
            conf_A2B = []
            conf_B2A = []
            for cr1, cr2 in zip(stats.credible_region( post1, self.conf), stats.credible_region( post2, self.conf)):
                if len(cr1) and len(cr2):
                    i, u = stats.geometric_overlap( cr1, cr2, nside, degrees=True )
                else:
                    i = u = 0
                conf_I.append( i )
                conf_U.append( u )

                conf_A2B.append( np.sum(post2[cr1]) )
                conf_B2A.append( np.sum(post1[cr2]) )

            ### things involving regions of fixed size
            area_I = []
            area_U = []
            area_A2B = []
            area_B2A = []
            for j in [int(np.ceil(1.0*area/pixarea)) for area in self.area]:
                i, u = stats.geometric_overlap( post1argsort[:j+1][:], post2argsort[:j+1][:], nside, degrees=True )
                area_I.append( i )
                area_U.append( u )

                area_A2B.append( np.sum(post2[post1argsort[:j+1]]) )
                area_B2A.append( np.sum(post1[post2argsort[:j+1]]) )

            self.comp["%s|%s"%(fits1,fits2)] = {'fidelity'  : stats.fidelity( post1, post2 ), ### fidelity
                                                'dTheta'    : np.arccos( stats.cos_dtheta( t1, p1, t2, p2) ), ### angular separation of mAP
                                                'conf'      : {'intersection' : conf_I, ### intersection and overlap of confidence regions
                                                               'union'        : conf_U,
                                                               'frAtoB'       : conf_A2B, ### confidence from post2 within CR defined by post1
                                                               'frBtoA'       : conf_B2A,
                                                              },
                                                'area'      : {'intersection' : area_I, ### intersection and overlap of first X deg2
                                                               'union'        : area_U,
                                                               'frAtoB'       : area_A2B, ### confidence from post2 within area defined by post1
                                                               'frBtoA'       : area_B2A,
                                                              },
                                               }

            ### plot on figures
            color = getColor.next()
            label = "%s - %s"%(self.texlabels[fits1], self.texlabels[fits2])

            # cr intersection, union, and ratio
            cri_ax.semilogy( self.conf, conf_I, color=color, label=label )
            cri_fig.fig.text(0.10+0.02, 0.93-0.05*ind, label, color=color, ha='left', va='top')

            cru_ax.semilogy( self.conf, conf_U, color=color, label=label )
            cru_fig.fig.text(0.10+0.02, 0.93-0.05*ind, label, color=color, ha='left', va='top')

            crr_ax.plot( self.conf, np.array(conf_I)/np.array(conf_U), color=color, label=label )
            crr_fig.fig.text(0.10+0.02, 0.93-0.05*ind, label, color=color, ha='left', va='top')

            # area intersection, union, and ratio
            ari_ax.loglog( self.area, area_I, color=color, label=label )
            ari_fig.fig.text(0.10+0.02, 0.93-0.05*ind, label, color=color, ha='left', va='top')

            aru_ax.loglog( self.area, area_U, color=color, label=label )
            aru_fig.fig.text(0.10+0.02, 0.93-0.05*ind, label, color=color, ha='left', va='top')

            arr_ax.semilogx( self.area, np.array(area_I)/np.array(area_U), color=color, label=label )
            arr_fig.fig.text(0.10+0.02, 0.93-0.05*ind, label, color=color, ha='left', va='top')

            # contained probabilities
            crc_ax.plot( self.conf, conf_A2B, color=color, linestyle='-', label=label )
            crc_tx.plot( self.conf, conf_B2A, color=color, linestyle='--', label=label )
            crc_fig.fig.text(0.50, 0.93-0.05*(ind+1), label, color=color, ha='center', va='top')

            arc_ax.semilogx( self.area, area_A2B, color=color, linestyle='-', label=label )
            arc_tx.semilogx( self.area, area_B2A, color=color, linestyle='--', label=label )
            arc_fig.fig.text(0.50, 0.93-0.05*(ind+1), label, color=color, ha='center', va='top')

        ### save json
        jsonname = "%s_compStats%s.js"%(self.label, self.tag)
        if verbose:
            print "  "+jsonname
        self.jsComp = Json( self.comp,
                          self.output_dir,
                          self.output_url,
                          graceid    = self.graceid,
                          graceDbURL = self.graceDbURL,
                          upload     = self.upload,
                        ).saveAndUpload( jsonname )

        ### save figures

        # cr intersection
        cri_ax.set_xlim(xmin=0.0, xmax=1.0)
        cri_ax.set_ylim(ymin=0.0)
        cri_ax.set_xlabel('confidence')
        cri_ax.set_ylabel('intersection [deg$^2$]')

        figname = "%s_CRIntersection%s.%s"%(self.label, self.tag, self.figtype)
        if verbose:
            print "  "+figname
        self.comp['cri'] = cri_fig.saveAndUpload( figname )

        plt.close( cri_fig.fig )

        # cr union
        cru_ax.set_xlim(xmin=0.0, xmax=1.0)
        cru_ax.set_ylim(ymin=0.0)
        cru_ax.set_xlabel('confidence')
        cru_ax.set_ylabel('union [deg$^2$]')

        figname = "%s_CRUnion%s.%s"%(self.label, self.tag, self.figtype)
        if verbose:
            print "  "+figname
        self.comp['cru'] = cru_fig.saveAndUpload( figname )

        plt.close( cru_fig.fig )

        # cr ratio
        crr_ax.set_xlim(xmin=0.0, xmax=1.0)
        crr_ax.set_ylim(ymin=0.0)
        crr_ax.set_xlabel('confidence')
        crr_ax.set_ylabel('intersection / union')

        figname = "%s_CRRatio%s.%s"%(self.label, self.tag, self.figtype)
        if verbose:
            print "  "+figname
        self.comp['crr'] = crr_fig.saveAndUpload( figname )

        plt.close( crr_fig.fig )

        # cr contained
        crc_ax.set_xlim(xmin=0.0, xmax=1.0)
        crc_ax.set_ylim(ymin=0.0, ymax=1.0)
        crc_ax.set_xlabel('$\int\limits_{\mathrm{CR}_A} d\Omega\, p_A(\Omega|\mathrm{data})$')
        crc_ax.set_ylabel('$\int\limits_{\mathrm{CR}_A} d\Omega\, p_B(\Omega|\mathrm{data})$')

        crc_tx.set_xlim(crc_ax.get_xlim())
        crc_tx.set_ylim(crc_ax.get_ylim())
        crc_tx.set_xlabel('$\int\limits_{\mathrm{CR}_B} d\Omega\, p_B(\Omega|\mathrm{data})$', visible=True)
        crc_tx.set_ylabel('$\int\limits_{\mathrm{CR}_B} d\Omega\, p_A(\Omega|\mathrm{data})$', visible=True)

        crc_tx.xaxis.tick_top()
        crc_tx.xaxis.set_label_position('top')
        crc_tx.yaxis.tick_right()
        crc_tx.yaxis.set_label_position('right')

        figname = "%s_CRContained%s.%s"%(self.label, self.tag, self.figtype)
        if verbose:
            print "  "+figname
        self.comp['crc'] = crc_fig.saveAndUpload( figname )

        plt.close( crc_fig.fig )

        # area intersection
        ari_ax.set_xlim(xmin=self.area[0], xmax=self.area[-1])
        ari_ax.set_ylim(ymin=0.0)
        ari_ax.set_xlabel('area [deg$^2$]')
        ari_ax.set_ylabel('intersection [deg$^2$]')

        figname = "%s_AreaIntersection%s.%s"%(self.label, self.tag, self.figtype)
        if verbose:
            print "  "+figname
        self.comp['ari'] = ari_fig.saveAndUpload( figname )

        plt.close( ari_fig.fig )

        # area union
        aru_ax.set_xlim(xmin=self.area[0], xmax=self.area[-1])
        aru_ax.set_ylim(ymin=0.0)
        aru_ax.set_xlabel('area [deg$^2$]')
        aru_ax.set_ylabel('union [deg$^2$]')

        figname = "%s_AreaUnion%s.%s"%(self.label, self.tag, self.figtype)
        if verbose:
            print "  "+figname
        self.comp['aru'] = aru_fig.saveAndUpload( figname )

        plt.close( aru_fig.fig )

        # area ratio
        arr_ax.set_xlim(xmin=self.area[0], xmax=self.area[-1])
        arr_ax.set_ylim(ymin=0.0)
        arr_ax.set_xlabel('area [deg$^2$]')
        arr_ax.set_ylabel('intersection / union')

        figname = "%s_AreaRatio%s.%s"%(self.label, self.tag, self.figtype)
        if verbose:
            print "  "+figname
        self.comp['arr'] = arr_fig.saveAndUpload( figname )

        plt.close( arr_fig.fig )

        # area contained
        arc_ax.set_xlim(xmin=self.area[0], xmax=self.area[-1])
        arc_ax.set_ylim(ymin=0.0, ymax=1.0)
        arc_ax.set_xlabel('$\int\limits_{\mathrm{CR}_A} d\Omega$')
        arc_ax.set_ylabel('$\int\limits_{\mathrm{CR}_A} d\Omega\, p_B(\Omega|\mathrm{data})$')

        arc_tx.set_xlim(arc_ax.get_xlim())
        arc_tx.set_ylim(arc_ax.get_ylim())
        arc_tx.set_xlabel('$\int\limits_{\mathrm{CR}_B} d\Omega$', visible=True)
        arc_tx.set_ylabel('$\int\limits_{\mathrm{CR}_B} d\Omega\, p_A(\Omega|\mathrm{data})$', visible=True)

        arc_tx.xaxis.tick_top()
        arc_tx.xaxis.set_label_position('top')
        arc_tx.yaxis.tick_right()
        arc_tx.yaxis.set_label_position('right')

        figname = "%s_AreaContained%s.%s"%(self.label, self.tag, self.figtype)
        if verbose:
            print "  "+figname
        self.comp['arc'] = arc_fig.saveAndUpload( figname )

        plt.close( arc_fig.fig )

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

    def __str__(self):
        '''
        generate html document as string
        '''
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
        head.link( href        = "https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css",
                   rel         = "stylesheet",
                   integrity   = "sha384-BVYiiSIFeK1dGmJRAkycuHAHRg32OmUcww7on3RYdg4Va+PmSTsz/K68vbdEjh4u",
                   crossorigin = "anonymous",
                 )

        ### other header information
        head.meta(name="description", content="a comparison of %s"%(', '.join(self.fitsnames)))
        head.meta(name="author", content=getpass.getuser()) ### whoever ran this is the author

        head.title(self.label)

        #----------------
        ### build body
        body = htmldoc.body()

        ### add navbar
        div1 = body.nav(klass='navbar navar-inverse navebar-fixed-top').div(klass='container').div(klass='navbar-header') ### first div, has links to snglFITS pages
        for fitsname in self.fitsnames:
            label = self.labels[fitsname]
            div1.a(label, klass='navbar-brand', href=os.path.join(self.output_url, "%s%s.html"%(label, self.tag)))

        body.hr

        div1 = body.nav(klass='navbar navar-inverse navebar-fixed-top').div(klass='container').div(klass='navbar-header') ### second div, has links to sections

        #---
        ### add sections
        sections = body.div(klass='container')

        ### top level summary
        if hasattr(self, 'fitsdata'): ### we must have called self.readFITS(), make_json, make_cumulative_json
            sections.hr
            div = sections.div(id='summary', klass='container')
            div.h1('Summary', id='summary')
            div1.a('Summary', klass='navbar-brand', href='#summary') ### add to top-level navigation bar

            for fitsname in self.fitsnames:
                row = div.div(klass='row') ### row for fitsname, nside, entropy, and information

                fitsbasename = os.path.basename(fitsname)

                row.div(klass='col-md-3').a(fitsbasename, href=os.path.join(self.output_url, fitsbasename))
                row.div(klass='col-md-2').p('nside = %d'%self.fitsdata[fitsname]['nside'])

        ### add mollweide plots
        if hasattr(self, 'mollweide'): ### we must have called make_mollweide
            sections.hr
            div = sections.div(id='mollweide', klass='container')
            div.h1('Mollweide Projections', id='mollweide')
            div1.a('Mollweide Projections', klass='navbar-brand', href='#mollweide') ### add to top-level navigation bar

            ### iterate over coordinates
            width = '575' ### FIXME: hard coded figure width isn't great...
            for coord, label in [("C","Equatorial"), ("E", "Geographic")]:
                div.h2(label+' coordinates', id="mollweide"+label)

                row = div.div(klass='row')
                row.div(klass='col-md-6').img(src=self.mollweide[coord], width=width)
                row.div(klass='col-md-6').img(src=self.mollweide[coord+' ann'], width=width)

        ### add confidence region plots
        if hasattr(self, 'CR'): ### must have called make_confidence_regions
            sections.hr
            div = sections.div(id='confidence regions', klass='container')
            div.h1('Confidence Regions', id='confidence')
            div1.a('Confidence Regions', klass='navbar-brand', href='#confidence') ### add to top-level navigation bar

            row = div.div(klass='row')

            ### put in the figures
            width = '550' ### FIXME: hard coding width isn't great...
            col = row.div(klass='col-md-6')
            col.img(src=self.CR['size'], width=width)
            col.img(src=self.CR['dTheta'], width=width)
            col.img(src=self.CR['modes'], width=width)

        ### add line-of-sight sanity checks and dT marginals
        if hasattr(self, 'dT') and hasattr(self, 'los'): ### must have called make_los and make_dT
            sections.hr
            div = sections.div(id='line of sight', klass='container')
            div.h1('Time Delay Marginals and Line-of-Sight Frame', id='timeDelay')
            div1.a('Time Delay', klass='navbar-brand', href='#timeDelay') ### add to top-level navigation bar

            for ifo1, ifo2 in self.ifo_pairs: ### add blurb for each IFO pair considered
                div.h2('%s - %s'%(ifo1, ifo2))
                ifos = "%s%s"%(ifo1, ifo2)

#                row = div.div(klass='row')
#                for fits1, fits2 in self.fits_pairs: ### formatting could be improved...
#                    col = row.div(klass='col-md-5')
#                    col.h3('%s - %s'%(self.labels[fits1], self.labels[fits2]))
#                    col.p().raw_text('&Delta;&theta;<sub>MAP</sub> = %.2f&deg;'%((self.dT[ifos][fits1]['thetaMAP']-self.dT[ifos][fits2]['thetaMAP'])*180/np.pi))
#                    col.p().raw_text('F(&Delta;t<sub>%s%s</sub>) = %.3f'%(ifo1,ifo2,self.dT[ifos]["%s|%s"%(fits1,fits2)]['fidelity']))

                row = div.div(klass='row')
                ### second col contains time-delay marginals
                col = row.div(klass='col-md-8')
                width = '700' ### FIXME: hard coding width isn't great...
#               col.img(src=self.dT[ifos]['fig'], width=width)
                col.img(src=self.dT[ifos]['ann fig'], width=width)
                width = '900' ### FIXME: hard coding width isn't great...
                col.img(src=self.los[ifos], width=width)

        ### add comparison statistics
        if hasattr(self, 'comp'): ### must have called make_comparison
            sections.hr
            div = sections.div(id='comparison statistics', klass='container')
            div.h1('Comparison Statistics', id='comparison')
            div1.a('Comparison Statistics', klass='navbar-brand', href='#comparison')

            row = div.div(klass='row')
            for fits1, fits2 in self.fits_pairs:
                col = row.div(klass='col-md-6')
                col.h3('%s - %s'%(self.labels[fits1], self.labels[fits2]))

                d = self.comp["%s|%s"%(fits1,fits2)]
                col.p().raw_text('F(&Omega;) = %.3f'%(d['fidelity']))
                col.p().raw_text('&Delta;&theta;<sub>MAP</sub> = %.2f&deg;'%(d['dTheta']*180/np.pi))

                if hasattr(self, 'dT'):
                    r = col.div(klass='row')
                    for ifo1, ifo2 in self.ifo_pairs:
                        c = r.div(klass='col-md-4')
                        c.h4('%s - %s'%(ifo1, ifo2))
                        ifos = "%s%s"%(ifo1, ifo2)

                        c.p().raw_text('&Delta;&theta;<sub>MAP</sub> = %.2f&deg;'%((self.dT[ifos][fits1]['thetaMAP']-self.dT[ifos][fits2]['thetaMAP'])*180/np.pi))
                        c.p().raw_text('F(&Delta;t<sub>%s%s</sub>) = %.3f'%(ifo1,ifo2,self.dT[ifos]["%s|%s"%(fits1,fits2)]['fidelity']))

            ### comparison statistics based on confidence regions
            row = div.div(klass='row')
            row.h2('comparisons based on confidence regions')

            col = row.div(klass='col-md-8')
            width = '550' ### FIXME: hard coding width isn't great...
            col.img(src=self.comp['cri'], width=width)
            col.img(src=self.comp['cru'], width=width)
            col.img(src=self.comp['crr'], width=width)
            col.img(src=self.comp['crc'], width=width)
           
            ### comparison statistics based on areas of fixed size
            row = div.div(klass='row')
            row.h2('comparisons based on areas of fixed size')

            col = row.div(klass='col-md-8')
            width = '550' ### FIXME: hard coding width isn't great...
            col.img(src=self.comp['ari'], width=width)
            col.img(src=self.comp['aru'], width=width)
            col.img(src=self.comp['arr'], width=width)
            col.img(src=self.comp['arc'], width=width)

        #----------------
        ### print document and return
        return str(doc)
