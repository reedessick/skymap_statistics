description = """generate some sanity checks based on basic triangulation"""
usage = "sanitycheck_maps.py [--options]"
author = "Reed Essick (reed.essick@ligo.org)"

import triangulate

from optparse import OptionParser

#=================================================

parser = OptionParser(usage=usage, description=description)

parser.add_option("-v", "--verbose", default=False, action="store_true")

opts, args = parser.parse_args()

#=================================================

"""
desired functionality:
    distributions of dtheta from known poles
    2D histograms rotated into the known-pole coordinate frame (mollweide projections?)
    "distance from bulk" measures? -> need to remember not only the line-of-sight but also the overhead directions for detectors

all of these can be accomplished by rotating to a specific coordinate system (depending on the goal) and then measuring a distribution. 
    -> rotate2pole
Either 1D or 2D histograms should cover this, with weighting by the probability distribution at that point.
Although, if you can set up some sort of interpolation algorithm, that would be wicked pissah.
"""



