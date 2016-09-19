# skymap_statistics
simple module and exectutables to quantify skymaps and compare them. Also contains plotting and visualization routines, as well as a script that builds friendly html documents summarizing this information (intended to be used with GraceDb).

--------------------------------------------------

# executables

All executables have descriptive help strings, but we briefly describe them here:

  - computation of statistics
    - analyze_maps.py 
        - computes statistics for each map separately, like confidence region sizes and skymap entropy. Can process multiple maps at a time, and simply iterates.
    - compare_maps.py 
        - computes statistics that compare maps in a pairwise manner. This includes the fidelity and the overlap of confidence regions. Can process an arbitrary number of maps, reporting each possible pair. 
    - antenna_averages.py 
        - computes statistics for the antenna patterns (like the values at the map's maximum a posteriori) from individual IFOs. Useful when considering the ratio of the expected strains in each IFO. Can also write the antenna patterns into a FITS file for latter use.

  - plotting and visualization
    - plot_maps.py 
        - provides basic mollweide projections as well as annotations based on triangulation and known detectors. Can process multiple maps at a time, generating a separate figure for each, as well as stacking contours of all maps together if requested.
    - overlay_maps.py :
        - generates pair-wise mollweide proejctions of contour overlays of all maps provided, producing a separate figure for each pair. Supports the same annotations as plot_maps.py.
    - plot_dt.py :
        - generates marginal distributions of the time-delay between two IFOs. Can process multiple maps and multiple IFO pairs at the same time. Generates a separate figure for each map and IFO pair provided, and will stack the posteriors for each IFO pair if requested. Particularly useful for determining whether the triangulation rings selected agree.
    - sanitycheck_maps.py : 
        - generates cartesisan projections of the maps provided in special frames defined by either the line-of-sight between detectors or the zenith of detectors. Produces a separate figure for each map provided and for each line-of-sight or zenith specified. Will also stack contours of posteriors if requested. Line-of-sight projections are particularly useful when determining if the support from different maps is similar along the ring.

  - coordinate transformations
    - rotate_maps.py
         - simply rotates a HEALPix representation from Equatorial to Geographic coordinates (or vice versa). 

  - summary html documents
    - snglFITShtml.py
        - builds an html document summarizing a single map. Will *not* process more than one map at a time. This is meant to provide a visual summary of most statistics within GraceDb (a link that references static URLs within GraceDb's file system). At this time, it is *not* recommended that the casual user attempt to run this script.

# libraries and dependencies

This module requires LALSuite (computes the Greenwich Mean Sidereel Time), Healpy, and Numpy. All other necessary code should be included within this repo, and most modules should be easy to use if users wish to create advanced analyses or plots not provided by the standard tools. 

--------------------------------------------------

# installation

There is no formal packaging for this module, at least not yet. However, we do provide a basic bash script that will set up your paths. 

To get the code, do

> cd ${srcdir}
> git clone https://github.com/reedessick/skymap_statistics.git

To set up paths (this will need to be done each time you want to run the library)

> cd ${srcdir}/skymap_statistics
> . setup.sh

You can confirm your set-up with

> which plot_maps.py

and

> plot_maps.py -h

--------------------------------------------------

# examples

Here, we provide an example of a complicated mollweide projection in Geographic coordinates, including annotations for the LHO antenna patterns and triangulation rings, assuming we have a single map originally specified in Equatorial coordinates. This should give users a sense of how to interact with some of the routines.

Assume we have a FITS file called EquatorialMap.fits.gz, which is in Equatorial coordinates and corresponds to an event at 1126259462.391. First, we need to rotate the FITS file into Geographic coordinates before we plot it.

> gps=1126259462.391
> rotate_map.py -v -S C -T E ${gps} EquatorialMap.fits.gz GeographicMap.fits.gz

This should produce a file called GeographicMap.fits.gz. Now, we need to compute the antenna patterns (as a FITS file)

> antenna_averages.py -o H -o L -c E --gps ${gps} --tag Geographic --write-fits GeographicMap.fits.gz

This should produce two files: H-antennaNorm_Geographic.fits and L-antennaNorm_Geographic.fits. It will also print some statistics about the antenna patterns to stdout. Because it will be useful down the road, I also normalize the antenna patterns so the map sums to one.

> python -c "import healpy as hp ; import numpy as np ; post=hp.read_map('H-antennaNorm_Geographic.fits') ; hp.write_map('H-antennaNorm_Geographic.fits', post/np.sum(post)) ; post=hp.read_map('L-antennaNorm_Geographic.fits') ; hp.write_map('L-antennaNorm_Geographic.fits', post/np.sum(post))"

Now, we want to actually generate the mollweide projection. There are a lot of moving parts here, so please reference plot_maps.py's help string for details on what each option actually controls. Note, we'll generate a few plots here with a few different options.

Let's start by just creating a mollweide projection in Equatorial coordinates.

> plot_maps.py -v myMap,EquatorialMap.fits.gz --projection "astro mollweide" --coord C --tag Equatorial --gps ${gps}

which produces myMap_Equatorial.png. We can create the analogous projection in Geographic coordinates (with continents outlined)

> plot_maps.py -v myMap,GeographicMap.fits.gz --projection "mollweide" --coord E --tag Geographic --gps ${gps} --continents

which is saved to myMap_Geographic.png. Now, let's make a fancier Geographic projection with annotations for a few interesting directions

> plot_maps.py -v myMap,GeographicMap.fits.gz --projection "mollweide" --coord E --tag Geographic --gps ${gps} --continents --line-of-sight HL --line-of-sight HV --line-of-sight LV --zenith H --zenith L --zenith V

Let's repeat that, but add a marker for a special spot (assume spherical coordinates of Theta=160deg and Phi=-5deg) and the associated triangulation rings

> plot_maps.py -v myMap,GeographicMap.fits.gz --projection "mollweide" --coord E --tag Geographic --gps ${gps} --continents --line-of-sight HL --line-of-sight HV --line-of-sight LV --zenith H --zenith L --zenith V --time-delay HL --time-delay HV --time-delay LV --time-delay-Dec-RA 160 -6 --time-delay-degrees --marker-Dec-RA 160 -6 --marker-degrees

This produces a pretty full plot, which shows all the basic triangulation information. For most maps, this should suffice to demonstrate the structure and patterns inherent within the posterior distribution.

Now I show how to make a mollweide projection which has the antenna pattern contours (we normalized the antenna pattern map to get the contours to work well).

> plot_maps.py -v LHO,H-antennaNorm_Geographic.fits LLO,L-antennaNorm_Geographic.fits --projection "mollweide" --coord E --tag Geographic --gps ${gps} --continents --stack-posteriors --stack-posteriors-background GeographicMap.fits.gz --stack-posteriors-levels 0.10 --stack-posteriors-levels 0.25 --stack-posteriors-levels 0.50 --stack-posteriors-levels 0.75 --stack-posteriors-levels 0.90 --line-of-sight HL --line-of-sight HV --line-of-sight LV --zenith H --zenith L --zenith V

Note, we use GeographicMap.fits.gz as the "background" in this plot so it is still shown as a heatmap. We then use H-antennaNorm_Geographic.fits as the map, and the --stack-posteriors option makes this a set of contours plotted on top of the heatmap in the figure called stackedPosterior_Geographic.png.

If we wish to produce the equivalent projection in Equatorial coordinates, we can with

> antenna_averages.py -o H -o L -c C --gps ${gps} --tag Equatorial --write-fits EquatorialMap.fits.gz

> python -c "import healpy as hp ; import numpy as np ; post=hp.read_map('H-antennaNorm_Equatorial.fits') ; hp.write_map('H-antennaNorm_Equatorial.fits', post/np.sum(post)) ; post=hp.read_map('L-antennaNorm_Equatorial.fits') ; hp.write_map('L-antennaNorm_Equatorial.fits', post/np.sum(post))"

> plot_maps.py -v LHO,H-antennaNorm_Equatorial.fits LLO,L-antennaNorm_Equatorial.fits --projection "astro mollweide" --coord C --tag Equatorial --gps ${gps} --stack-posteriors --stack-posteriors-background EquatorialMap.fits.gz --stack-posteriors-levels 0.10 --stack-posteriors-levels 0.25 --stack-posteriors-levels 0.50 --stack-posteriors-levels 0.75 --stack-posteriors-levels 0.90 --line-of-sight HL --line-of-sight HV --line-of-sight LV --zenith H --zenith L --zenith V

We'll also quickly generate a marginal distribution for the delay-time between LHO and LLO. We can do this with either map; notice how the output is identical.

> plot_dt.py -v myMap,GeographicMap.fits.gz --time-delay HL --coord E --gps ${gps} --tag Geographic

> plot_dt.py -v myMap,EquatorialMap.fits.gz --time-delay HL --coord C --gps ${gps} --tag Equatorial
