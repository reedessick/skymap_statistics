#!/usr/bin/python

description = """the author is far too well educated to spend his time on data entry, so this automates the comparisons he put together of various skymaps. Assumes that all relevant information is attached to GraceDb (as it should be)."""
usage = """autosummary.py [--options] gracedid graceid graceid ..."""
author = "Reed Essick (reed.essick@ligo.org)"

import os
import sys
import subprocess as sp

import numpy as np
import stats
import triangulate

from ligo.gracedb.rest import GraceDb

from optparse import OptionParser

#=================================================

def plot_cmd( fitsfile, radec=None ):
    """
    return the command string for plotting.
    if radec!=None:
        ra, dec = radec
    """
    figname = "%s/%s.png"%(os.path.dirname(fitsfile), os.path.basename(fitsfile).split(".")[0])
    string = "bayestar_plot_allsky %s -o %s"%(fitsfile, figname)
    if radec:
        bayestar_plot_allsky += " --radec %.5f %.5f"%radec
    return string, figname

###

def analyze_cmd( fitsfile, radec=None ):
    """
    return the command string for analyzing
    """
    string = "python analyze_maps.py -v -d -H -c 0.10 -c 0.25 -c 0.50 -c 0.75 -c 0.90 -c 0.95 -c 0.99 %s,%s"%(os.path.basename(fitsfile).split(".")[0], fitsfile)
    if radec:
        ra, dec = radec
        phi = ra*180/np.pi
        theta = (dec - 0.5*np.pi)*180/np.pi
        string += " --pvalue %.5f,%.5f --searched-area %.5f,%.5f"(theta, phi, theta, phi)
    return string

###

def sanitycheck_cmd( fitsfile, geocent, outdir=".", los="H,L" ):
    """
    return the command string for sanity checking
    """
    string = "python sanitycheck_maps.py -v -L %s -c C -T %.6f -o %s -t %s -m -p -g %s,%s"%(los, geocent, outdir, os.path.basename(fitsfile).split(".")[0], os.path.basename(fitsfile).split(".")[0], fitsfile)
    return string

###

def sanityoverlay_cmd( fitsfiles, geocent, outdir=".", los="H,L" ):
    """
    return the command string for sanity overlay plots
    """
    string = "python sanitycheck_maps.py -v -L %s -c C -T %.6f -o %s -p -g -C %s"%(los, geocent, outdir, " ".join( "%s,%s"%(os.path.basename(f).split(".")[0], f) for f in fitsfiles ) )
    return string

###

def compare_cmd( fitsfile, fitsFile ):
    """
    return the command string for comparing
    """
    string = "python compare_maps.py -v -d --dMAP --fidelity --structural-similarity -s 0.10 -s 0.25 -s 0.50 -s 0.75 -s 0.90 -s 0.95 -s 0.99 -c 0.10 -c 0.25 -c 0.50 -c 0.75 -c 0.90 -c 0.95 -c 0.99 %s,%s %s,%s"%(os.path.basename(fitsfile).split(".")[0], fitsfile, os.path.basename(fitsFile).split(".")[0], fitsFile)
    return string

###

def overlay_cmd( fitsfile, fitsFile, outdir="." ):
    """
    return the command string for overlaying
    """
    string = "python overlay_maps.py -v -c 0.10 -c 0.50 -c 0.90 -c 0.99 -o %s %s,%s %s,%s"%(outdir, os.path.basename(fitsfile).split(".")[0], fitsfile, os.path.basename(fitsFile).split(".")[0], fitsFile)
    return string

###

def compile_cmd( docname ):
    """
    return the command string for latex compilation
    """
    string = "pdflatex %s"%(docname)
    return string

#=================================================

def data2latex( event, fitsfiles ):
    """
    returns a latex document as a string
    This is the work-horse of the script, where all the tex formatting and annoying crap like that is applied.
    """
    graceid = event['graceid']
    ### start off string with a preamble, etc
    string = r"""\documentclass{article}

\usepackage{fullpage}
\usepackage{graphicx}
\usepackage{hyperref}

\title{
skymap comparison for \href{https://gracedb.ligo.org/events/view/%s}{%s}
}
\author{
autosummary.py\footnote{\url{https://github.com/reedessick/skymap_statistics}}
}

\begin{document}

\maketitle
"""%(graceid, graceid)

    ### section describing top-level attributes
    group = event['group']
    pipeline = event['pipeline']
    if event.has_key('search'):
        search = event['search']
    else:
        search = ''

    gpstime = float(event['gpstime'])
    created = event['created']

    far = float(event['far'])

    instruments = event['instruments']
    labels = ",".join(event['labels'].keys())

    string += r"""
\section{event description}

\begin{center}
    \begin{tabular}{ccc| c c | c | c | c}
        group & pipeline & search & gpstime & created & far [Hz] & instruments & labels \\
        \hline
        %s    & %s       & %s     & %.9f    & %s      & %.5e     & %s          & %s 
    \end{tabular}
\end{center}
"""%(group, pipeline, search, gpstime, created, far, instruments, labels)

    fitsorder = sorted( key for key in fitsfiles.keys() if "fit" in key ) ### gets rid of sanity overlay stuff

    ### section describing each map
    string += r"""
\section{descriptions of individual maps}
"""

    for fitsname in fitsorder:
        plot = fitsfiles[fitsname]['plot']
        string += r"""
\subsection{%s}

\begin{center}
    \includegraphics[width=1.0\textwidth]{%s}
\end{center}
"""%(fitsname.replace("_", "\_"), os.path.basename(plot))

        ### analyze
        analyze = fitsfiles[fitsname]['analyze']
        string += r"""
%s
"""%analyze2string( analyze )

        ### sanitycheck
        sanity = fitsfiles[fitsname]['sanitycheck']
        string += r"""
%s
"""%sanity2string( sanity )

    ### section comparing map
    if len(fitsorder) > 1:
        string += r"""
\section{comparison of maps}
"""
        ### HOW DO I BUILD A TABLE NICELY FROM THE EXISTING DATA? 
        ### need to set up an array and fill it in, iterating through the data
        ### then, and only then, should I attempt to print a string-formatted-thing to the document
 
        ### compare
        string += r"""
%s
"""%(compare2string( fitsorder, fitsfiles ))

        ### overlay
        string += r"""
%s
"""%(overlay2string( fitsorder, fitsfiles ))

        ### sanityoverlay plot
        string += r"""
%s
"""%sanityoverlay2string( fitsfiles['sanityoverlay'] )

    ### finish document
    string += r"""
\end{document}"""

    return string

###

def analyze2string( analyze ):
    CR = {}
    data = {}
    for line in analyze.split("\n"):
        if "CR" in line:
            fields = line.strip().split()
            cr = float(fields[0])*1e-2
            if "=" in line:
                if CR.has_key(cr):
                    CR[cr][fields[3]] = float(fields[-2])
                else:
                    CR[cr] = {fields[3]:float(fields[-2])}
            else: ### disjoint regions
                if CR.has_key(cr):
                    CR[cr]['disjoint regions'] = [float(l) for l in line.split("(")[-1].split(")")[0].split(", ")]
                else:
                    CR[cr] = {'disjoint regions':[float(l) for l in line.split("(")[-1].split(")")[0].split(", ")]}
        elif "=" in line:
            fields = line.strip().split("=")
            data[fields[0].strip()] = float(fields[1].split()[0])
           
    ### this could use improvement! 
    ### hard code this based on conditionals for what is present
    string = r"""
\begin{itemize}"""
    if data.has_key("nside"):
        string += r"""
    \item{nside : %d}"""%(data['nside'])
    if data.has_key("pixarea"):
        string += r"""
    \item{pixarea : %.5f deg}"""%(data['pixarea'])
    if data.has_key("entropy"):
        string += r"""
    \item{entropy : %.5f deg$^2$}"""%(data['entropy'])

    ### need to include searched area and/or pvalue?

    string += r"""
\end{itemize}
"""

    string += r"""
\begin{center}
    \begin{tabular}{c | c c c}
        confidence & size [deg$^2$] & disjoint regions [deg$^2$] & max\{$\delta \theta$\} [deg] \\
        \hline"""
    for cr in sorted(CR.keys()):
        size = CR[cr]['size']
        disjoint_regions = CR[cr]['disjoint regions']
        dtheta = CR[cr]['max(dtheta)']
        string += r"""
        %.1f %s & %.5f & %s & %.5f \\"""%(round(cr*1e2, 1), "\%", size, ", ".join("%.5f"%l for l in disjoint_regions), dtheta)
    string = string[:-2] ### remove the trailing "\\"
    string += r"""
    \end{tabular}
\end{center}
"""

    return string

###

def sanity2string( sanity ):
    string = ""

    for line in sanity.split("\n"):
        if ("los" in line) and (".png" in line):
            figname = os.path.basename(line.strip())
            break
    else:
        raise ValueError("could not find a figure file...")

    string += r"""
\begin{center}
    \includegraphics[width=1.0\textwidth]{%s}
\end{center}
"""%(figname)

    miD = []
    for line in sanity.split("\n"):
        if "mutualinformationDistance" in line:
           miD.append( float(line.strip().split()[-1]) )
    if len(miD)!=2:
        raise ValueError("could not find mutualInformationDistance...")

    string += r"""
\begin{center}
    \begin{tabular}{c | c c}
        & unrotated & rotated \\
        \hline
        mutual information distance $\equiv \frac{I(\theta,\phi)}{H(\theta,\phi)}$ & %.6f & %.6f
    \end{tabular}
\end{center}
"""%tuple(miD)

    return string

###

def sanityoverlay2string( sanityoverlay ):
    for line in sanityoverlay.split("\n"):
        if ("los" in line) and (".png" in line):
            figname = os.path.basename(line.strip())
            break
    else:
        raise ValueError("could not find a figure file...")

    string = r"""
\begin{center}
    \includegraphics[width=1.0\textwidth]{%s}
\end{center}
"""%(figname)

    return string

###

def compare2string( fitsorder, fitsfiles ):

    string = ""

    string += r"""
\begin{center}
    \begin{tabular}{c| %s }
        $\delta \theta_{MAP}$ & %s \\
        \hline"""%((len(fitsorder)-1)*"c", " & ".join(l.replace("_","\_") for l in fitsorder[1:]))
    for fitsname in fitsorder[:-1]:
        string += r"""
        %s  & \\"""%(fitsname.replace("_","\_"))
    string = string[:-2]
    string += r"""
    \end{tabular}
\end{center}
"""

    string += r"""
\begin{center}
    \begin{tabular}{c| %s }
        fidelity & %s \\
        \hline"""%((len(fitsorder)-1)*"c", " & ".join(l.replace("_","\_") for l in fitsorder[1:]))
    for fitsname in fitsorder[:-1]:
        string += r"""
        %s & \\"""%(fitsname.replace("_","\_"))
    string = string[:-2]
    string += r"""
    \end{tabular}
\end{center}
"""

    string += r"""
\begin{center}
    \begin{tabular}{c| %s }
        structural similarity & %s \\
        \hline"""%((len(fitsorder)-1)*"c", " & ".join(l.replace("_","\_") for l in fitsorder[1:]))
    for fitsname in fitsorder[:-1]:
        string += r"""
        %s & \\"""%(fitsname.replace("_","\_"))
    string = string[:-2]
    string += r"""
    \end{tabular}
\end{center}
"""

    for cr in [0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99]:
        string += r"""
\begin{center}
    \begin{tabular}{c| %s }
        confidence region : %.1f %s intersection / union & %s \\
        \hline"""%((len(fitsorder)-1)*"c", cr*1e2, "\%", " & ".join(l.replace("_","\_") for l in fitsorder[1:]))
        for fitsname in fitsorder[:-1]:
            string += r"""
        %s & \\"""%(fitsname.replace("_","\_"))
        string = string[:-2]
        string += r"""
    \end{tabular}
\end{center}
"""
        string += r"""
\begin{center}
    \begin{tabular}{c| %s }
        spot check : %.1f %s & %s \\
        \hline"""%((len(fitsorder))*"c", cr*1e2, "\%"," & ".join(l.replace("_","\_") for l in fitsorder))
        for fitsname in fitsorder:
            string += r"""
        %s & \\"""%(fitsname.replace("_","\_"))
        string = string[:-2]
        string += r"""
    \end{tabular}
\end{center}
"""

    return string

### 

def overlay2string( fitsorder, fitsfiles ):
    string = r"""
\begin{center}
    \begin{tabular}{c| %s }
        overlay & %s \\
        \hline"""%((len(fitsorder)-1)*"c", " & ".join(l.replace("_","\_") for l in fitsorder[1:]))
    for fitsname in fitsorder[:-1]:
        string += r"""
%s \\"""%(fitsname.replace("_","\_"))
    string = string[:-2]
    string += r"""
    \end{tabular}
\end{center}
"""

    return string

#=================================================

parser = OptionParser(description=description, usage=usage)

parser.add_option("-v", "--verbose", default=False, action="store_true")

parser.add_option('-G', '--gracedb_url', default=None, type="string")

parser.add_option("-o", "--output-dir", default=".", type="string")

parser.add_option("-c", "--compile", default=False, action="store_true")

opts, args = parser.parse_args()

if not len(args):
    raise ValueError("please supply at least one graceid as an input argument")

#=================================================

if opts.gracedb_url:
    if opts.verbose:
        print( "conecting to GraceDb : %s"%(opts.gracedb_url) )
    gracedb = GraceDb( opts.gracedb_url )
else:
    if opts.verbose:
        print("connecting to GraceDb")
    gracedb = GraceDb()

#=================================================

### iterate through each event, creating as separate document for each
for graceid in args:
    if opts.verbose:
        print( "processing : %s"%(graceid) )

    ### set up directory
    outdir = "%s/%s"%(opts.output_dir, graceid)
    if not os.path.exists( outdir ):
        os.makedirs( outdir )

    ### get event data
    if opts.verbose:
        print( "retrieving event data" )
    event = gracedb.event( graceid ).json()
    geocent = float(event['gpstime'])
    instruments = event['instruments']

    if instruments.replace("1","")!="H,L":
        raise ValueError("don't know how to processes events with instruments=%s"%(instruments))

    ### pull down fits files
    if opts.verbose:
        print( "pulling down data for" )
    files = sorted( gracedb.files( graceid ).json().keys() )
    fitsfiles = dict( (filename,{}) for filename in files if filename.strip(".gz").endswith(".fits") )
    for fits in fitsfiles.keys():
        filename = "%s/%s"%(outdir, fits)
        if opts.verbose:
            print( "\t%s"%(filename) )
        file_obj = open(filename, "w")
        file_obj.write( gracedb.files( graceid, fits ).read() )
        file_obj.close()
        fitsfiles[fits]['path'] = filename
 
    ### plot fits files
    if opts.verbose:
        print( "plotting skymaps" )
    for fitsfile in fitsfiles.keys():
        path = fitsfiles[fitsfile]['path']
        cmd, outfilename = plot_cmd( path, radec=None )
        if opts.verbose:
            print( "\t%s"%cmd )
        sp.Popen( cmd.split() ).wait()
        fitsfiles[fitsfile]['plot'] = outfilename

    ### analyse maps
    if opts.verbose:
        print( "analyzing maps" )
    for fitsfile in fitsfiles.keys():
        path = fitsfiles[fitsfile]['path']
        cmd = analyze_cmd( path, radec=None )
        if opts.verbose:
            print( "\t%s"%cmd )
        result = sp.Popen( cmd.split(), stdout=sp.PIPE ).communicate()[0]
        fitsfiles[fitsfile]['analyze'] = result

    files = fitsfiles.keys()

    ### sanity check maps
    if opts.verbose:
        print( "sanity checking maps" )
    for fitsfile in fitsfiles.keys():
        path = fitsfiles[fitsfile]['path']
        cmd = sanitycheck_cmd( path, geocent, outdir=outdir, los="H,L" ) ### IFO THING MAY NEED UPDATING...
        if opts.verbose:
            print( "\t%s"%cmd )
        result = sp.Popen( cmd.split(), stdout=sp.PIPE ).communicate()[0]
        fitsfiles[fitsfile]['sanitycheck'] = result

    if opts.verbose:
        print( "sanity overlay" )
    cmd = sanityoverlay_cmd( [fitsfiles[fitsfile]['path'] for fitsfile in files], geocent, outdir=outdir, los="H,L" )
    if opts.verbose:
        print( "\t%s"%cmd )
    result = sp.Popen( cmd.split(), stdout=sp.PIPE ).communicate()[0]
    fitsfiles['sanityoverlay'] = result

    ### compare maps
    if opts.verbose:
        print( "comparing maps" )
    for ind, fitsfile in enumerate( files ):
        path1 = fitsfiles[fitsfile]['path']

        for fitsFile in files[ind+1:]:
            path2 = fitsfiles[fitsFile]['path']

            cmd = compare_cmd( path1, path2 )
            if opts.verbose:
                print( "\t%s"%cmd )
            result = sp.Popen( cmd.split(), stdout=sp.PIPE ).communicate()[0]
            fitsfiles[fitsfile]['compare : %s'%fitsFile] = result
            fitsfiles[fitsFile]['compare : %s'%fitsfile] = result

    ### overlay maps
    if opts.verbose:
        print( "overlay maps" )
    for ind, fitsfile in enumerate( files ):
        path1 = fitsfiles[fitsfile]['path']

        for fitsFile in files[ind+1:]:
            path2 = fitsfiles[fitsFile]['path']

            cmd = overlay_cmd( path1, path2, outdir=outdir )
            if opts.verbose:
                print( "\t%s"%cmd )
            result = sp.Popen( cmd.split(), stdout=sp.PIPE ).communicate()[0]
            fitsfiles[fitsfile]['overlay : %s'%fitsFile] = result
            fitsfiles[fitsFile]['overlay : %s'%fitsfile] = result

    ### write latex document
    docname = "%s/summary.tex"%(outdir)
    if opts.verbose:
        print( "writing : %s"%(docname) )
    doc = open(docname, "w")
    doc.write( data2latex( event, fitsfiles ) )
    doc.close()

    if opts.compile:
        cmd = compile_cmd( docname )
        if opts.verbose:
            print( "\t%s"%cmd )
        sp.Popen(cmd.split()).wait()



