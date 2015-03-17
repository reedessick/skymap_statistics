usage="""a module to compute basic statistics of skymaps"""

import healpy as hp
import numpy as np


#=================================================
#
# general helper methods
#
#=================================================
def cos_dtheta(theta1, phi1, theta2, phi2):
	"""
	computes the angular separation between two points
	support arrays assuming they all have the same shape
	"""
	return np.cos(theta1)*np.cos(theta2) + np.sin(theta1)*np.sin(theta2)*np.cos(phi1-phi2)

#=================================================
#
# methods involving manipulating a posterior
#
#=================================================
def rankmap(posterior, npix=None):
        """
        converts a posterior into a rank map. Small ranks correspond to large posterior weight
        WARNING: rankmaps do not sum to unity!
        """
	if not npix:
	        npix = len(posterior)
        rankmap = np.empty(npix,int)
        rankmap[posterior.argsort()] = np.arange(npix)[::-1]
        return rankmap

###
def resample(posterior, new_nside, nest=False):
	"""
	creates a new posterior with len=hp.nside2npix(new_nside)
		if new_nside > nside: we assign each new pixel an equal fraction of the parent pixel's value
		if new_nside < nside: we assign each new pixel the sum of the contained pixels
	"""
	return hp.ud_grade(posterior, new_nside, power=-2)

###
def __to_cumulative(posterior):
	"""
	returns a map corresponding to cumulative probabilities at each pixel
	assumes ``greedy binning'' algorithm
	"""
	cum = np.zeros_like( posterior )
	c = 0.0
	for i in posterior.argsort()[::-1]:
		c += posterior[i]
		cum[i] = c
	return cum

###
def credible_region(posterior, conf):
	"""
	returns a list of pixels that correspond the minimum credible region
	"""
	if isinstance(conf, (int,float)):
		conf = np.array([conf])
	elif not isinstance(conf, np.ndarray):
		conf = np.array(conf)
	if np.any(conf < 0) or np.any(conf > 1):
		raise ValueError("conf must be between 0 and 1")
	
	cum = __to_cumulative(posterior)
	
	ind = np.arange(len(posterior))
	return [ind[cum<=c] for c in conf] ### return indecies corresponding to confidence levels

###
def p_value(posterior, theta, phi, nside=None):
	"""
	computes the p-value at which a given point was found
	"""
	if not nside:
		nside = hp.npix2nside(len(posterior))
	return np.sum(posterior[posterior>=posterior[hp.ang2pix(nside, theta, phi)]])

###
def estang(posterior, nside=None):
	"""
	returns the position associated with the maximum of the posterior
	"""
	if not nside:
		nside = hp.npix2nside(len(posterior))
	return hp.pix2ang(nside, posterior.argmax())

###
def searched_area(posterior, theta, phi, nside=None, nest=False, degrees=False):
	"""
	computes the searched area given a location
	"""
	if not nside:
		nside = hp.npix2nside(len(posterior))
	ipix = hp.ang2pix(nside, theta, phi, nest=nest)
	return np.sum(posterior>=posterior[ipix])*hp.nside2pixarea(nside, degrees=degrees)

###
def est_cos_dtheta(posterior, theta, phi):
	"""
	returns the angular separation between the maximum of the posterior and theta, phi
	"""
	t, p = estang(posterior)
	return cos_dtheta(theta, phi, t, p)

###
def min_cos_dtheta(posterior, theta, phi, nside=None, nest=False):
	"""
	computes the maximum angular separation between any point in the area with p > p(theta, phi) and the estimated position
	"""
	if not nside:
		nside = hp.npix2nside(len(posterior))
	ipix = hp.ang2pix(nside, theta, phi, nest=nest)
	## get all ang included in this area
	thetas, phis = hp.pix2ang(nside, np.arange(len(posterior))[posterior >= posterior[ipix]])

	t, p = estang(posterior)

	return np.min(cos_dtheta(thetas, phis, t, p))

###
def min_all_cos_dtheta(pix, nside, nest=False):
	"""
	computes the maximum angular separation between any two pixels within pix=[ipix,ipix,...]
	"""
	min_c = 1.0
	thetas, phis = hp.pix2ang(nside, pix)
	for i, (t1, p1) in enumerate(zip(thetas, phis)[:-1]):
		t2 = thetas[i+1:] 
		p2 = phis[i+1:]
		c = np.amin(cos_dtheta(t1,p1,t2,p2))
		if c < min_c:
			min_c = c
	return min_c


###
def num_modes(posterior, theta, phi, nside=None, nest=False):
	"""
	computes the number of modes in the area bounded by theta, phi
	"""
	npix = len(posterior)
	if not nside:
		nside = hp.npix2nside(npix)
	pix = list(np.arange(npix)[posterior>=posterior[hp.ang2pix(nside, theta, phi, nest=nest)]]) ## get list of pixels
	return len( __into_modes(nside, pix, nest=nest) )

### 
def size_modes(posterior, theta, phi, nside=None, nest=False, degrees=False):
	"""
	computes the sizes of modes in the area bounded by theta, phi
	"""
	npix = len(posterior)
	if not nside:
		nside = hp.npix2nside(npix)
	pix = list(np.arange(npix)[posterior>=posterior[hp.ang2pix(nside, theta, phi, nest=nest)]])
	return [len(_)*hp.nside2pixarea(nside, degrees=degrees) for _ in __into_modes(nside, pix)]

def __into_modes(nside, pix, nest=False):
	"""
	divides the list of pixels (pix) into simply connected modes
	"""
	### establish an array representing the included pixels
	npix = hp.nside2npix(nside)

	truth = np.zeros((npix,),bool)
	truth[pix] = True

	pixnums = np.arange(npix) ### convenient array we establish once

	modes = []
	while truth.any():
		ipix = pixnums[truth][0] ### take the first pixel
		truth[ipix] = False ### remove it from the global set
		mode = [ipix]
		to_check = [ipix] ### add it to the list of things to check

		while len(to_check): # there are pixels in this mode we have to check
			ipix = to_check.pop() # take one pixel from those to be checked.

			for neighbour in hp.get_all_neighbours(nside, ipix, nest=nest):# get neighbors as rtheta, rphi

				if neighbour == -1: ### when neighbour == -1, there is no corresponding pixel in this direction
					pass
				# try to find pixel in skymap
				elif truth[neighbour]: ### pixel in the set and has not been visited before
					truth[neighbour] = False ### remove neighbour from global set
#					truth[neighbour] = 0 ### remove neighbour from global set
					mode.append( neighbour ) ### add to this mode
					to_check.append( neighbour ) ### add to list of things to check
				else: ### pixel not in the set or has been visited before
					pass 
		modes.append( mode )

	return modes

###
def entropy(posterior, base=2.0):
	"""
	computes the shannon entropy in the posterior
	we compute the entropy with base=base
		base=2 => bits
		base=e => nats
	"""
	posterior = posterior[posterior > 0.0] ### log(0)*0 -> 0
	return -np.sum( np.log(posterior)*posterior)/np.log(base)

###
def information(posterior, base=2.0):
	"""
	computes the shannon entropy in the posterior
	we compute the entropy with base=base
		base=2 => bits
		base=e => nats
	"""
	n = len(posterior)
	posterior = posterior[posterior > 0.0] ### log(0)*0 -> 0
	logbase = np.log(base)
	return (np.log(n) + np.sum( np.log(posterior)*posterior) )/np.log(base)

#=================================================
#
# methods for comparing two skymaps
#
#=================================================
def mse(posterior1, posterior2):
	"""
	computes the mean square error between the two posteriors
		sum (p1 - p2)**2
	"""
	return 1.0*np.sum( (posterior1-posterior2)**2 )/len(posterior1)

###
def peak_snr(posterior1, posterior2):
	"""
	computes the peak signal-to-noise ratio between the two posteriors
		max(posterior1)/mse , max(posterior2)/mse
	"""
	_mse = mse(posterior1, posterior2)**0.5
	return np.max(posterior1)/_mse, np.max(posterior2)/_mse

###
def fidelity(posterior1, posterior2):
	"""
	computes the fidelity between the two posteriors
		sum (p1*p2)**0.5
	"""
	return np.sum( (posterior1*posterior2)**0.5 )

###
def KLdivergence(posterior1, posterior2, base=2.0):
	"""
	computes the Kullback-Leibler divergence
		sum log(p1/p2) * p1
	"""
	truth = posterior1>0
	return np.sum( np.log(posterior1[truth]/posterior2[truth])*posterior1[truth] )/np.log(base)

###
def symmetric_KLdivergence(posterior1, posterior2, base=2.0):
	"""
	computes the symmetric Kullback-Leibler divergence
		sum log(p1/p2)*(p1 - p2)
	"""
	return KLdivergence(posterior1, posterior2, base=base) + KLdivergence(posterior2, posterior2, base=base)
#	return np.sum( np.log(posterior1/posterior2)*(posterior1 - posterior2) )/np.log(base)

###
def structural_similarity(posterior1, posterior2, c1=(0.01)**2, c2=(0.03)**2):
	"""
	computes the structural similarity
		(2*m1*m2+c1)*(2*v12 + c2) / (m1**2 + m2**2 + c1)*(v1 + v2 + c2)
	dynamic range of our pixels is 1
	m1 => mean(p1)
	v1 => var(p1)
	v12=> covar(p1,p2)
	"""
	m1 = np.mean(posterior1)
	m2 = np.mean(posterior2)
	covar = np.cov(posterior1, posterior2)
	return ( (2*m1*m2+c1)*(2*covar[0,1] + c2) ) / ( (m1**2 + m2**2 + c1)*(covar[0,0] + covar[1,1] + c2) )

### 
def pearson(posterior1, posterior2):
	"""
	computes the pearson correlation coefficient
	"""
	covar = np.cov(posterior1, posterior2)
	return covar[0,1]/(covar[0,0]*covar[1,1])**0.5

###
def dot(posterior1, posterior2):
	"""
	takes the "inner product" of the two posteriors and returns the cos(theta) between them
	"""
	return np.sum(posterior1*posterior2)/(np.sum(posterior1**2)*np.sum(posterior2**2))**0.5

###
def geometric_overlap(pix1, pix2, nside, degrees=False):
        """
        computes the amount of area in the intersection and union of confidence regions from p1 and p2 defined by p_value
        """
        npix = hp.nside2npix(nside)
        pixarea = hp.nside2pixarea(nside, degrees=degrees)

	posterior1 = np.zeros((npix,),int)
	posterior2 = np.zeros((npix,),int)

	posterior1[pix1] = 1
	posterior2[pix2] = 1

        intersection = np.sum( posterior1*posterior2 )
        return intersection*pixarea, (np.sum(posterior1+posterior2) - intersection)*pixarea

