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

#def plot_cmd( fitsfile, radec=None ):
#    """
#    return the command string for plotting.
#    if radec!=None:
#        ra, dec = radec
#    """
#    figname = "%s/%s.png"%(os.path.dirname(fitsfile), os.path.basename(fitsfile).split(".")[0])
#    string = "bayestar_plot_allsky %s -o %s"%(fitsfile, figname)
#    if radec:
#        string += " --radec %.5f %.5f"%radec
#    return string, figname

def plot_cmd( fitsfile, label, outdir=".", radec=None, color_map="Reds" ):
    string = "python plot_maps.py -v -o %s --color-map %s %s,%s"%(outdir, color_map, label, fitsfile)
    return string, "%s/%s.png"%(outdir, label)

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

def sanitycheck_cmd( fitsfile, geocent, outdir=".", los="H,L", color_map="Reds" ):
    """
    return the command string for sanity checking
    """
    label = os.path.basename(fitsfile).split(".")[0]
    string = "python sanitycheck_maps.py -v -L %s -c C -T %.6f -o %s -t %s -m -p -g --color-map %s %s,%s"%(los, geocent, outdir, label, color_map, label, fitsfile)
    return string

###

def sanityoverlay_cmd( fitsfiles, geocent, outdir=".", los="H,L" ):
    """
    return the command string for sanity overlay plots
    """
    string = "python sanitycheck_maps.py -v -L %s -c C -T %.6f -o %s -p -g -C %s"%(los, geocent, outdir, " ".join( "%s:%s,%s"%(g,l, f) for g,f,l in fitsfiles ) )
    return string

###

def compare_cmd( (gid, fitsfile, label), (gID, fitsFile, laBel) ):
    """
    return the command string for comparing
    """
    string = "python compare_maps.py -v -d --dMAP --fidelity --structural-similarity --symKL -s 0.10 -s 0.25 -s 0.50 -s 0.75 -s 0.90 -s 0.95 -s 0.99 -c 0.10 -c 0.25 -c 0.50 -c 0.75 -c 0.90 -c 0.95 -c 0.99 %s:%s,%s %s:%s,%s"%(gid, label, fitsfile, gID, laBel, fitsFile)
    return string

###

def overlay_cmd( (gid, fitsfile, label), (gID, fitsFile, laBel), outdir="." ):
    """
    return the command string for overlaying
    """
    string = "python overlay_maps.py -v -c 0.10 -c 0.50 -c 0.90 -c 0.99 -o %s %s:%s,%s %s:%s,%s"%(outdir, gid, label, fitsfile, gID, laBel, fitsFile)
    return string

###

def compile_cmd( docname ):
    """
    return the command string for latex compilation
    """
    string = "pdflatex %s"%(docname)
    return string

#=================================================

def data2latex( event, fitsfiles, landscape=False, neighbors=[]):
    """
    returns a latex document as a string
    This is the work-horse of the script, where all the tex formatting and annoying crap like that is applied.
    """
    graceid = event['graceid']
    ### start off string with a preamble, etc
    if landscape:
        docclass="[landscape]{article}"
    else:
        docclass="{article}"
    string = r"""\documentclass%s

\usepackage{fullpage}
\usepackage{graphicx}
\usepackage{hyperref}
\usepackage{multirow}

\begin{document}

\title{
skymap comparison for \href{https://gracedb.ligo.org/events/view/%s}{%s}
}
\author{
autosummary.py\footnote{\url{https://github.com/reedessick/skymap_statistics}}
}

\maketitle

\abstract{
This document was generated and compiled automatically and has not been validated before publication.
}

"""%(docclass, graceid, graceid)

    ### section describing top-level attributes
    group = event['group']
    pipeline = event['pipeline']
    if event.has_key('search'):
        search = event['search']
    else:
        search = ''

    gpstime = float(event['gpstime'])
    created = event['created']

    try:
        far = float(event['far'])
    except:
        far = np.nan

    instruments = event['instruments']
    labels = ",".join(event['labels'].keys())

    string += r"""
\section{event description}

\begin{center}
    \begin{tabular}{c | ccc| c c | c | c | c}
        GraceID & group & pipeline & search & gpstime & created & far [Hz] & instruments & labels \\
        \hline
        %s      & %s    & %s       & %s     & %.9f    & %s      & %.5e     & %s          & %s \\"""%(graceid, group, pipeline, search, gpstime, created, far, instruments, labels)

    for e in neighbors:
        egraceid = e['graceid']
        egroup = e['group']
        epipeline = e['pipeline']
        if e.has_key('search'):
            esearch = e['search']
        else:
            search = ''

        egpstime = float(e['gpstime'])
        ecreated = e['created']

        try:
            efar = float(e['far'])
        except:
            efar = np.nan

        einstruments = e['instruments']
        elabels = ",".join(e['labels'].keys())

        string += r"""
        %s      & %s    & %s       & %s     & %.9f    & %s      & %.5e     & %s          & %s \\"""%(egraceid, egroup, epipeline, esearch, egpstime, ecreated, efar, einstruments, elabels)

    string = string[:-2]
    string += r"""
    \end{tabular}
\end{center}
"""
    fitsorder = [ key for key in fitsfiles.keys() if isinstance(key, tuple) and ("fit" in key[1]) ] ### gets rid of sanity overlay stuff
    fitsorder.sort( key=lambda l: l[1] )
    fitsorder.sort( key=lambda l: l[0] )

    if not len(fitsorder):
        string += r"""
\end{document}"""
        return string

    string += r"""
\begin{itemize}"""
    for gid, fitsname in fitsorder:
        string += r"""
    \item{%s : %s}"""%(gid, fitsfiles[(gid,fitsname)]['label'].replace("_","\_"))
    string += r"""
\end{itemize}
"""

    ### section describing each map
    string += r"""
\section{descriptions of individual maps}
"""

    for fitsname in fitsorder:
        label = fitsfiles[fitsname]['label']
        plot = fitsfiles[fitsname]['plot']
        string += r"""
\subsection{%s:%s}

\begin{center}
    \includegraphics[width=0.45\textwidth]{%s}
    \includegraphics[width=0.45\textwidth]{%s}
\end{center}
"""%(fitsname[0], label.replace("_", "\_"), "../%s/%s"%tuple(plot.split("/")[-2:]), sanity2plot( fitsfiles[fitsname]['sanitycheck'] ))

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

    ### glossary
    string += r"""
\newpage
\section{glossary}

\subsection{descritpion of individual maps}

\subsubsection{mollweide projcetion}

Standard Mollweide projection of our skymap in Celestial coordinates.

\subsubsection{line-of-sight projection with marginals}

A cartesian projection of the skymap in a frame defined by the line-of-sight between LHO and LLO.
Because most of the localization comes from triangulation, we see that the distributions are approximately
separable in this frame.
In fact, we can measure the time-of-arrival difference between LHO and LLO based on the marignal distribution over $\theta$.

\subsubsection{nside, pixarea}

Basic parameters of the HEALPix decomposition used for this particular map.

\subsubsection{entropy}

A measure of the uncertainty associated with a map based on the Shannon entropy.
We normalize this to a fraction of the sky (measured in deg$^2$) to avoid any ambiguity due to different levels of pixelization.

\begin{equation}
    H(p) = A_\mathrm{pix} e^{-\sum_\mathrm{pix} p \ln p}
\end{equation}

\subsubsection{Confidence Regions size, No. disjoint regions, $\mathrm{max}\{\delta\theta\}$}

We define the confidence region $R_c$ cooresponding to confidence $c$ as the minimum area that contains at a probability greater or equal to $c$

\begin{equation}
    R_c := \min\limits_{\Sigma} \int\limits_\Sigma\mathrm{d}\Omega\ \left|\ \int_\Sigma\mathrm{d}\Omega\, p \geq c \right.
\end{equation}

Typically, this is computed via a pixelated skymap.
We include pixels in order of decreasing probability until we reach a cumulative probability $\geq c$.
In this way, we uniquely identify a set of pixels associated with $R_c$.
The size is just the corresponding solid angle associated with this set.
The No. disjoint regions is the number of separate regions, which do not touch, contained in the map.
Typically, for 2-detector networks, we have 2 disjoint regions (one overhead and one underneath the detectors).
$\mathrm{max}\{\delta\theta\}$ is the largest angular separation between any two pixels contained in the set.

\subsection{Mutual Information}

We define the mutual information between $\theta$ and $\phi$ as the Kullback-Leibler divergence from the joint distribution to the product of the marginals.

\begin{equation}
    I(\theta,\phi) = D_\mathrm{KL}(p(\theta,\phi)||p(\theta)p(\phi)) = \int\sin\theta\mathrm{d}\theta\mathrm{d}\phi\, p(\theta,\phi) \ln \frac{p(\theta,\phi)}{p(\theta)p(\phi)}
\end{equation}

where $p(\theta) = \int\mathrm{d}\phi\, p(\theta,\phi)$ and $p(\phi) = \int\sin\theta\mathrm{d}\theta\, p(\theta,\phi)$. 
Note, this measure represents the amount of information lost by approximating the true joint distribution by the product of the marginals.
It is \textit{not} invarient under general coordinate transformations, in particular under rotations. 
We use this to our advantage and compute $I$ for the original map and the map rotated into the line-of-sight frame defined by the two detectors.
If our intuition from triangulation is correct, then $I$ in the rotated frame should be significantly smaller (typically a factor of $\sim$10) than $I$ in Celestial coordinates.

\subsection{comparison of maps}

\subsubsection{$\delta\theta_\mathrm{MAP}$}

Simply the angular separation between the maximum a posteriori points of two maps. 

\subsubsection{Fidelity}

A measure of the distance between two probability distributions, the fidelity is defined as

\begin{equation}
    F(p, q) = \int\mathrm{d}\Omega \sqrt{p \cdot q} \in [0, 1]
\end{equation}

with $F=1$ iff $p(\Omega)=q(\Omega)\ \forall\Omega$.
Essentially, it is a modified dot product that tells us how similar to distributions are.

\subsubsection{Structural Similarity}

The structural similarity attempts to model how the human eye distinguishes between images.

\begin{equation}
    \mathrm{SSI} = \left( \frac{2\mu_p\mu_q + c_1}{\mu_p^2 + \mu_q^2 + c_2}\right)\left( \frac{2\Sigma_{pq} + c_2}{\Sigma_{pp} + \Sigma_{qq} + c_2}\right) \in [-1, 1]
\end{equation}

where $\mu_{p} = \sum_\mathrm{pix} p / N_\mathrm{pix}$ and $\Sigma_{pq} = \sum_\mathrm{pix} (p-\mu_p)(q-\mu_q) / N_\mathrm{pix}$, with anagogous definitions for the other terms.
$c_1$ and $c_2$ are regularizing constants that prevent the statistic from diverging if the averages and/or variances vanish. 
These are taken to be much smaller than the expected range of $\mu$ and $\Sigma$ and therefore have little effect on the statistic in practice.

\subsubsection{Symmetrick Kullback-Leibler divergence}

The Kullback-Leibler divergence is a non-symmetric distance measure between to distributions.
It can be thought of as the information lost when approximating the distribution $p$ with the distribution $q$.

\begin{equation}
    D_\mathrm{KL}(p||q) = \int\mathrm{d}\Omega\, p \ln \frac{p}{q} \in [0,\infty)
\end{equation}

We symmeterize this because we do not have an a priori preference between maps, which results in the symmetric KL divergence

\begin{equation}
    D_\mathrm{symKL}(p,q) = D_\mathrm{KL}(p||q) + D_\mathrm{KL}(q||p) \in [0,\infty)
\end{equation}

WARNING: this statistic is not particularly ``stable'' to slight changes in the pixelization and can often be infinite.

\subsubsection{Confidence Regions intersection, union, spotcheck}

We define confidence regions for each map separately, which defines a two sets of pixels.
The intersection and union of these sets of pixels are straightforward.

The \textit{spotcheck} is a non-symmetric measure of the probability contined in on map within a region defined by another map.
Precisely, it is 

\begin{equation}
    D_\mathrm{s}(p||q;c) = \int\limits_{R_c(q)}\mathrm{d}\Omega p\ \left|\ R_c(q) := \min\limits_{\Sigma} \int\limits_\Sigma\mathrm{d}\Omega\ \right|\ \int_\Sigma\mathrm{d}\Omega p \geq c 
\end{equation}

In more practical terms, if I observe all of the confidence region $R_c$ defined by map $q$, $D_\mathrm{s}$ defines the amount of probability contained in my observations based on map $p$.

\subsection{Contour overlays}

We provide overlays of the contours from separate maps in both Celestial and line-of-sight coordinates.
The line-of-sight distributions come with marginal distributions as well.
This helps identifiy whether the maps pick up the same triangulation rings and have the same wight around the rings.

"""

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

def sanity2plot( sanity ):
    string = ""
    for line in sanity.split("\n"):
        if ("los" in line) and (".png" in line):
            return "../%s/%s"%tuple(line.strip().split("/")[-2:])
    else:
        print sanity
        raise ValueError("could not find a figure file...")

###

def sanity2string( sanity ):
    string = ""

    miD = []
    for line in sanity.split("\n"):
        if "mutualinformationDistance" in line:
           miD.append( float(line.strip().split()[-1]) )
    if len(miD)!=2:
        print sanity
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
            figname = line.strip()
            break
    else:
        print sanityoverlay
        raise ValueError("could not find a figure file...")

    string = r"""
\begin{center}
    \includegraphics[width=0.9\textwidth]{../%s/%s}
\end{center}
"""%tuple(figname.split("/")[-2:])

    return string

###

def compare2spotcheck( compare, label1, label2, cr ):
    ### label1 = row -> integrated
    ### label2 = column -> defines CR
    look = "spotcheck %.3f %s"%(cr*100, "%")
    order = None
    for line in compare.split("\n"):
        if (label1 in line) and (label2 in line):
            order = tuple( line.split(" vs ") )
        if look in line:
            vals = line.split(":")[-1].strip(" (").strip(")").split(", ")
            break
    else:
        print compare
        raise ValueError("could not find spotcheck...")
    if not order:
        print compare
        raise ValueError("could not find order...")

    if order == (label1, label2):
        return vals[1] 
    else:
        return vals[0]

###

def compare2cr( compare, cr ):
    look = "%.3f %s CR"%(cr*100, "%")
    dat = []
    for line in compare.split("\n"):
        if look in line:
            dat.append( line.split(":")[-1].strip() )
    if len(dat) != 4:
        raise ValueError("something's wrong with dat's length...")
    return "%s$^\circ$ / %s$^\circ$"%(dat[-2].split()[-2], dat[-1].split()[-2])

###

def compare2structuralsimilarity( compare ):
    for line in compare.split("\n"):
        if "structural similarity" in line:
            return line.strip().split()[-1]
    else:
        print compare
        raise ValueError("could not find structural similarity...")

###

def compare2dthetaMAP( compare ):
    for line in compare.split("\n"):
        if "dtheta_MAP" in line:
            return line.strip().split()[-2]
    else:
        print compare
        raise ValueError("could not find dtheta_MAP...")

###

def compare2fidelity( compare ):
    for line in compare.split("\n"):
        if "fidelity" in line:
            return line.strip().split()[-1]
    else:
        print compare
        raise ValueError("could not find fidelity...")

###

def compare2symKL( compare ):
    for line in compare.split("\n"):
        if "symmetric KL divergence" in line:
            value = line.strip().split()[-1]
            if value == "inf":
                value = "$\infty$"
            return value
    else:
        print compare
        raise ValueError("could not find symmetric KL divergence...")

###

def compare2string( fitsorder, fitsfiles ):

    string = ""

    string += r"""
\subsection{$\delta \theta_{MAP}$}

\begin{center}
    \begin{tabular}{c| %s }
        $\delta \theta_{MAP}$ & %s \\
        \hline"""%((len(fitsorder)-1)*"c", " & ".join("%s:%s"%(l[0], fitsfiles[l]['label'].replace("_","\_")) for l in fitsorder[1:]))

    for ind, fitsname in enumerate(fitsorder[:-1]):
        label1 = fitsfiles[fitsname]['label']
        plusequal = r"""
        %s:%s"""%(fitsname[0],label1.replace("_","\_"))
        for i in xrange(ind):
            plusequal += r""" & -"""
        for fitsName in fitsorder[ind+1:]:
            label2 = fitsfiles[fitsName]['label']
            plusequal += r""" & %s$^\circ$"""%compare2dthetaMAP( fitsfiles[fitsname]['compare : %s_%s'%fitsName] )
        string += r"%s \\"%(plusequal)
    string = string[:-2]
    string += r"""
    \end{tabular}
\end{center}
"""

    string += r"""
\subsection{fidelity}

\begin{center}
    \begin{tabular}{c| %s }
        fidelity & %s \\
        \hline"""%((len(fitsorder)-1)*"c", " & ".join("%s:%s"%(l[0], fitsfiles[l]['label'].replace("_","\_")) for l in fitsorder[1:]))

    for ind, fitsname in enumerate(fitsorder[:-1]):
        label1 = fitsfiles[fitsname]['label']
        plusequal = r"""
        %s:%s"""%(fitsname[0], label1.replace("_","\_"))
        for i in xrange(ind):
            plusequal += r""" & -"""
        for fitsName in fitsorder[ind+1:]:
            label2 = fitsfiles[fitsName]['label']
            plusequal += r""" & %s"""%compare2fidelity( fitsfiles[fitsname]['compare : %s_%s'%fitsName] )
        string += r"%s \\"%(plusequal)

    string = string[:-2]
    string += r"""
    \end{tabular}
\end{center}
"""

    string += r"""
\subsection{structural similarity}

\begin{center}
    \begin{tabular}{c| %s }
        structural similarity & %s \\
        \hline"""%((len(fitsorder)-1)*"c", " & ".join("%s:%s"%(l[0], fitsfiles[l]['label'].replace("_","\_")) for l in fitsorder[1:]))

    for ind, fitsname in enumerate(fitsorder[:-1]):
        label1 = fitsfiles[fitsname]['label']
        plusequal = r"""
        %s:%s"""%(fitsname[0], label1.replace("_","\_"))
        for i in xrange(ind):
            plusequal += r""" & -"""
        for fitsName in fitsorder[ind+1:]:
            label2 = fitsfiles[fitsName]['label']
            plusequal += r""" & %s"""%compare2structuralsimilarity( fitsfiles[fitsname]['compare : %s_%s'%fitsName] )
        string += r"%s \\"%(plusequal)

    string = string[:-2]
    string += r"""
    \end{tabular}
\end{center}
"""

    string += r"""
\subsection{symmetric KL divergence}

\begin{center}
    \begin{tabular}{c| %s }
        symmetric KL divergence & %s \\
        \hline"""%((len(fitsorder)-1)*"c", " & ".join("%s:%s"%(l[0], fitsfiles[l]['label'].replace("_","\_")) for l in fitsorder[1:]))

    for ind, fitsname in enumerate(fitsorder[:-1]):
        label1 = fitsfiles[fitsname]['label']
        plusequal = r"""
        %s:%s"""%(fitsname[0], label1.replace("_","\_"))
        for i in xrange(ind):
            plusequal += r""" & -"""
        for fitsName in fitsorder[ind+1:]:
            label2 = fitsfiles[fitsName]['label']
            plusequal += r""" & %s"""%compare2symKL( fitsfiles[fitsname]['compare : %s_%s'%fitsName] )
        string += r"%s \\"%(plusequal)

    string = string[:-2]
    string += r"""
    \end{tabular}
\end{center}
"""

    string += r"""
\subsection{Confidence levels}"""

    for cr in [0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99]:
        string += r"""
\subsubsection{%.1f%s}"""%(cr*100, "\%")

        string += r"""
\begin{center}
    \begin{tabular}{c| %s }
        confidence region : %.1f%s intersection / union & %s \\
        \hline"""%((len(fitsorder)-1)*"c", cr*1e2, "\%", " & ".join("%s:%s"%(l[0], fitsfiles[l]['label'].replace("_","\_")) for l in fitsorder[1:]))

        for ind, fitsname in enumerate(fitsorder[:-1]):
            label1 = fitsfiles[fitsname]['label']
            plusequal = r"""
        %s:%s"""%(fitsname[0], label1.replace("_","\_"))
            for i in xrange(ind):
                plusequal += r""" & -"""
            for fitsName in fitsorder[ind+1:]:
                label2 = fitsfiles[fitsName]['label']
                plusequal += r""" & %s"""%compare2cr( fitsfiles[fitsname]['compare : %s_%s'%fitsName], cr )
            string += r"%s \\"%(plusequal)

        string = string[:-2]
        string += r"""
    \end{tabular}
\end{center}
"""

        string += r"""
\begin{center}
    \begin{tabular}{c| %s }
                             & \multicolumn{%d}{|c}{map defining confidence region} \\
        spot check : %.1f%s & %s \\
        \hline"""%((len(fitsorder))*"c", len(fitsorder), cr*1e2, "\%"," & ".join("%s:%s"%(l[0], fitsfiles[l]['label'].replace("_","\_")) for l in fitsorder))
 
        for ind, fitsname in enumerate(fitsorder):
            label1 = fitsfiles[fitsname]['label']
            plusequal = r"""
        %s:%s"""%(fitsname[0], label1.replace("_","\_"))
            for fitsName in fitsorder[:ind]:
                label2 = fitsfiles[fitsName]['label']
                plusequal += r""" & %s"""%compare2spotcheck( fitsfiles[fitsname]['compare : %s_%s'%fitsName], label1, label2, cr )
            plusequal += r""" & -"""
            for fitsName in fitsorder[ind+1:]:
                label2 = fitsfiles[fitsName]['label']
                plusequal += r""" & %s"""%compare2spotcheck( fitsfiles[fitsname]['compare : %s_%s'%fitsName], label1, label2, cr )
            string += r"%s \\"%(plusequal)

        string = string[:-2]
        string += r"""
    \end{tabular}
\end{center}
"""
    return string

### 

def overlay2plot( overlay ):
    for line in overlay.split("\n"):
        if ".png" in line:
            line = line.strip()
            return "../%s/%s"%tuple(line.strip().split("/")[-2:])
    else:
        print overlay
        raise ValueError("could not find a figure file...")

###

def overlay2string( fitsorder, fitsfiles ):
    string = r"""
\begin{center}
    \begin{tabular}{c| %s }
        overlay & %s \\
        \hline"""%((len(fitsorder)-1)*"c", " & ".join("%s:%s"%(l[0], fitsfiles[l]['label'].replace("_","\_")) for l in fitsorder[1:]))

    width = 0.9/len(fitsorder)
    for ind, fitsname in enumerate(fitsorder[:-1]):
        label1 = fitsfiles[fitsname]['label']
        plusequal = r"""
        %s:%s"""%(fitsname[0], label1.replace("_","\_"))
        for i in xrange(ind):
            plusequal += r""" &"""
        for fitsName in fitsorder[ind+1:]:
            label2 = fitsfiles[fitsName]['label']
            plusequal += r""" & \includegraphics[width=%.6f\textwidth]{%s}"""%(width, overlay2plot(fitsfiles[fitsname]['overlay : %s_%s'%fitsName]))
        string += r"%s \\"%(plusequal)

    string = string[:-2]
    string += r"""
    \end{tabular}
\end{center}
"""
    return string

#=================================================

parser = OptionParser(description=description, usage=usage)

parser.add_option("-v", "--verbose", default=False, action="store_true")

parser.add_option("-G", '--gracedb_url', default=None, type="string")

parser.add_option("-o", "--output-dir", default=".", type="string")

parser.add_option("-l", "--landscape", default=False, action="store_true", help="make latex doc lanscape")
parser.add_option("-c", "--compile", default=False, action="store_true", help="compile latex doc")

parser.add_option("-F", "--force", default=False, action="store_true", help="download all files and re-compute all statistics even if they already exist")

parser.add_option("-a", "--annotate-gracedb", default=False, action="store_true", help="upload a pdf to gracedb")

parser.add_option("-w", "--neighbor-window", default=None, type="float", help="search for neighbors within +/- neighbors_window and include any maps from those events in the comparison" )

parser.add_option("", "--neighbors-not-my-group", default=False, action="store_true", help="restrict neighbors to events of different group")
parser.add_option("", "--neighbors-not-my-pipeline", default=False, action="store_true", help="restrict neighbors to events from different pipeline")
parser.add_option("", "--neighbors-not-my-search", default=False, action="store_true", help="restrict neighbors to events from different search. Does nothing if the event has no search defined.")

parser.add_option("", "--color-map", default="cylon", type="string")

opts, args = parser.parse_args()

if not len(args):
    raise ValueError("please supply at least one graceid as an input argument")

notforce = not opts.force ### used often, so we compute it once

opts.compile = opts.compile or opts.annotate_gracedb

if opts.compile: ### needed for compilation
    cwd = os.getcwd()

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
        print( "\tretrieving event data" )
    event = gracedb.event( graceid ).json()
    geocent = float(event['gpstime'])
    instruments = event['instruments']

    if instruments.replace("1","")!="H,L":
        raise ValueError("don't know how to processes events with instruments=%s"%(instruments))

    ### pull down fits files
    files = sorted( gracedb.files( graceid ).json().keys() )
    fitsfiles = dict( ((graceid, filename),{}) for filename in files if filename.strip(".gz").endswith(".fits") )
    old = True ### records whether there is a new FITS
                ### used for sanityoverlay
    for _, fits in fitsfiles.keys():
        filename = "%s/%s"%(outdir, fits)
        if notforce and os.path.exists( filename ):
            if opts.verbose:
                print( "\t\t%s already exists"%(filename) )
        else:
            if opts.verbose:
                print( "\t\t%s"%(filename) )
            file_obj = open(filename, "w")
            file_obj.write( gracedb.files( graceid, fits ).read() )
            file_obj.close()
            old = False
        fitsfiles[(graceid, fits)]['path'] = filename
        fitsfiles[(graceid, fits)]['label'] = os.path.basename(filename).split(".")[0]
        fitsfiles[(graceid, fits)]['graceid'] = graceid

    if opts.neighbor_window > 0: ### perform neighbor search
        if opts.verbose:
            print( "\tsearching for neighbors within %.3f sec"%(opts.neighbor_window) )
        neighbors = [e for e in gracedb.events( "%.6f .. %.6f"%(geocent-opts.neighbor_window, geocent+opts.neighbor_window) ) if (e['graceid'][0] != 'H') and (e['graceid'] != graceid)]

        if opts.neighbors_not_my_group: ### only include events from different group
            neighbors = [e for e in neighbors if e['group'] != event['group']]
        elif opts.neighbors_not_my_pipeline: ### only include events from different pipeline
            neighbors = [e for e in neighbors if e['pipeline'] != event['pipeline']]
        elif opts.neighbors_not_my_search and event.has_key('search'): ### only include events from different search
            neighbors = [e for e in neighbors if (not e.has_key('search')) or (e['search']==event['search'])]

        for e in neighbors:
            if opts.verbose:
                print( "\t\tfound : %s"%(e['graceid']) )
            noutdir = "%s/%s"%(opts.output_dir, e['graceid'])
            if not os.path.exists( noutdir ):
                os.makedirs( noutdir )
            files = sorted( gracedb.files( e['graceid'] ).json().keys() )
            for fits in [filename for filename in files if filename.strip(".gz").endswith(".fits") ]:
                fitsfiles.update( {(e['graceid'],fits):{}} )
                filename = "%s/%s"%(noutdir, fits)
                if notforce and os.path.exists( filename ):
                    if opts.verbose:
                        print( "\t\t\t%s already exists"%(filename) )
                else:
                    if opts.verbose:
                        print( "\t\t\t%s"%(filename) )
                    file_obj = open(filename, "w")
                    file_obj.write( gracedb.files( e['graceid'], fits ).read() )
                    file_obj.close()
                    old = False
                fitsfiles[(e['graceid'], fits)]['path'] = filename
                fitsfiles[(e['graceid'], fits)]['label'] = os.path.basename(filename).split(".")[0]
                fitsfiles[(e['graceid'], fits)]['graceid'] = e['graceid']
    else:
        neighbors = []

    ### plot fits files
    if opts.verbose:
        print( "\tplotting skymaps" )
    for gid, fitsfile in fitsfiles.keys():
        path = fitsfiles[(gid, fitsfile)]['path']
        label = fitsfiles[(gid, fitsfile)]['label']
        cmd, outfilename = plot_cmd( path, label, outdir="%s/%s"%(opts.output_dir, gid), color_map=opts.color_map )
        if notforce and os.path.exists( outfilename ):
            if opts.verbose:
                print "\t\t%s already exists"%(outfilename)
        else:
            out = "%s.out"%(outfilename)
            err = "%s.err"%(out[:-4])
            if opts.verbose:
                print( "\t\t%s > %s, %s"%(cmd, out, err) )
            out_obj = open(out, "w")
            err_obj = open(err, "w")
            sp.Popen( cmd.split(), stdout=out_obj, stderr=err_obj ).wait()
            out_obj.close()
            err_obj.close()
        fitsfiles[(gid, fitsfile)]['plot'] = outfilename

    ### analyse maps
    if opts.verbose:
        print( "\tanalyzing maps" )
    for gid, fitsfile in fitsfiles.keys():
        path = fitsfiles[(gid, fitsfile)]['path']
        out = "%s.analyze.out"%(path)
        if notforce and os.path.exists( out ):
            if opts.verbose:
                print( "\t\t%s already exists"%(out) )
        else:
            err = "%s.err"%(out[:-4])
            cmd = analyze_cmd( path, radec=None )
            if opts.verbose:
                print( "\t\t%s > %s, %s"%(cmd, out, err) )
            out_obj = open(out, "w")
            err_obj = open(err, "w")
            sp.Popen( cmd.split(), stdout=out_obj, stderr=err_obj ).wait()
            out_obj.close()
            err_obj.close()
        out_obj = open( out , "r" )
        fitsfiles[(gid, fitsfile)]['analyze'] = out_obj.read()
        out_obj.close()

    files = fitsfiles.keys()
    files.sort(key=lambda l: l[1])
    ### sanity check maps
    if opts.verbose:
        print( "\tsanity checking maps" )
    for gid, fitsfile in fitsfiles.keys():
        path = fitsfiles[(gid, fitsfile)]['path']
        out = "%s.sanitycheck.out"%(path)
        if notforce and os.path.exists( out ):
            if opts.verbose:
                print( "\t\t%s already exists"%(out) )
        else:
            cmd = sanitycheck_cmd( path, geocent, outdir="%s/%s"%(opts.output_dir, gid), los="H,L", color_map=opts.color_map ) ### IFO THING MAY NEED UPDATING...
            err = "%s.err"%(out[:-4])
            if opts.verbose:
                print( "\t\t%s > %s, %s"%(cmd, out, err) )
            out_obj = open(out, "w")
            err_obj = open(err, "w")
            sp.Popen( cmd.split(), stdout=out_obj, stderr=err_obj ).wait()
            out_obj.close()
            err_obj.close()
        out_obj = open( out, "r" )
        fitsfiles[(gid, fitsfile)]['sanitycheck'] = out_obj.read()
        out_obj.close()

    ### sanity overlay
    out = "%s/sanityoverlay.out"%(outdir)
    if old and notforce:
        if opts.verbose:
            print( "\tnothing new to sanity overlay" )
    else:
        if opts.verbose:
            print( "\tsanity overlay" )
        cmd = sanityoverlay_cmd( [(gid, fitsfiles[(gid,fitsfile)]['path'], fitsfiles[(gid,fitsfile)]['label']) for gid,fitsfile in files], geocent, outdir=outdir, los="H,L" )
        err = "%s.err"%(out[:-4])
        if opts.verbose:
            print( "\t\t%s > %s, %s"%(cmd, out, err) )
        out_obj = open(out, "w")
        err_obj = open(err, "w")
        sp.Popen( cmd.split(), stdout=out_obj, stderr=err_obj ).wait()
        out_obj.close()
        err_obj.close()
    out_obj = open( out, "r" )
    fitsfiles['sanityoverlay'] = out_obj.read()
    out_obj.close()

    ### compare maps
    if opts.verbose:
        print( "\tcomparing maps" )
    for ind, (gid, fitsfile) in enumerate( files ):
        path1 = fitsfiles[(gid, fitsfile)]['path']
        label1 = fitsfiles[(gid, fitsfile)]['label']

        for gID, fitsFile in files[ind+1:]:
            path2 = fitsfiles[(gID, fitsFile)]['path']
            label2 = fitsfiles[(gID, fitsFile)]['label']

            out = "%s/%s_%s-%s_%s.compare.out"%(outdir, gid, label1, gID, label2)
            if notforce and os.path.exists( out ):
                if opts.verbose:
                    print( "\t\t%s already exists"%(out) )
            else:
                cmd = compare_cmd( (gid, path1, label1), (gID, path2, label2) )
                err = "%s.err"%(out[:-4])
                if opts.verbose:
                    print( "\t\t%s > %s, %s"%(cmd, out, err) )
                out_obj = open(out, "w")    
                err_obj = open(err, "w")    
                sp.Popen( cmd.split(), stdout=out_obj, stderr=err_obj ).wait()
                out_obj.close()
                err_obj.close()
            out_obj = open( out, "r" )
            result = out_obj.read()
            out_obj.close()
            fitsfiles[(gid,fitsfile)]['compare : %s_%s'%(gID,fitsFile)] = result
            fitsfiles[(gID,fitsFile)]['compare : %s_%s'%(gid,fitsfile)] = result

    ### overlay maps
    if opts.verbose:
        print( "\toverlay maps" )
    for ind, (gid, fitsfile) in enumerate( files ):
        path1 = fitsfiles[(gid, fitsfile)]['path']
        label1 = fitsfiles[(gid, fitsfile)]['label']

        for gID, fitsFile in files[ind+1:]:
            path2 = fitsfiles[(gID, fitsFile)]['path']
            label2 = fitsfiles[(gID, fitsFile)]['label']

            out = "%s/%s_%s-%s_%s.overlay.out"%(outdir, gid, label1, gID, label2)
            if notforce and os.path.exists( out ):
                if opts.verbose:
                    print( "\t\t%s already exists"%(out) )
            else:
                cmd = overlay_cmd( (gid, path1, label1), (gID, path2, label2), outdir=outdir )
                err = "%s.err"%(out[:-4])
                if opts.verbose:
                    print( "\t\t%s > %s, %s"%(cmd, out, err) )
                out_obj = open(out, "w")
                err_obj = open(err, "w")
                sp.Popen( cmd.split(), stdout=out_obj, stderr=err_obj ).communicate()[0]
                out_obj.close()
                err_obj.close()
            out_obj = open( out, "r" )
            result = out_obj.read()
            out_obj.close()
            fitsfiles[(gid, fitsfile)]['overlay : %s_%s'%(gID, fitsFile)] = result
            fitsfiles[(gID, fitsFile)]['overlay : %s_%s'%(gid, fitsfile)] = result

    #====================

    ### write latex document
    docname = "%s/summary.tex"%(outdir)
    if opts.verbose:
        print( "\twriting : %s"%(docname) )
    doc = open(docname, "w")
    doc.write( data2latex( event, fitsfiles, landscape=opts.landscape, neighbors=neighbors ) )
    doc.close()

    if opts.compile:
        os.chdir( outdir )
        cmd = compile_cmd( os.path.basename(docname) )
        if opts.verbose:
            print( "\t%s"%cmd )
        sp.Popen(cmd.split()).wait()
        os.chdir( cwd )

    if opts.annotate_gracedb:
        if opts.verbose:
            print( "\tuploading to GraceDb" )
        labels = [fitsfiles[(gid, fitsname)]['label'] for gid, fitsname in files ]
        lenlabels = len(labels)        
        if lenlabels > 2:
            labels = "%s, and %s"%(", ".join(labels[:-1]), labels[-1])
        elif lenlabels > 1:
            labels = " and ".join(labels)
        else:
            labels = labels[0]
        message = "skymap comparison for %s"%(labels)
        if opts.verbose:
            print( "\t%s"%(message) )
        gracedb.writeLog( graceid, message, filename=docname, tagname="sky_loc" )


